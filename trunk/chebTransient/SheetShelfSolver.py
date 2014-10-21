import numpy
import Chebyshev
import numpy.linalg as linalg
  
from optparse import OptionParser

from pprint import pprint

import os.path
import string

import scipy.optimize

from scipy.integrate import odeint 

def bcRes(xg,solver):
  
  (b,bx) = solver.computeB(xg)
  H = b/(1-solver.delta)
  u = solver.options.a*xg/H
  ux = (solver.delta/(8*solver.epsilon)*H)**solver.n
  Hx = (solver.options.a - ux*H)/u
  res = -solver.gamma*u**(1/solver.n)-H*(Hx-bx)
  return res

def dy_dx(y, x, solver):
  H = y[0]
  u = solver.options.a*x/H
  
  (b,bx) = solver.computeB(x)
  tauB = -solver.gamma*numpy.abs(u)**(1/solver.n-1)*u
  Hx = tauB/H + bx
  return [Hx] 

class SheetShelfSolver:
  def __init__(self, options):
    if(not os.path.exists(options.folder)):
      os.mkdir(options.folder)
  
    self.options = options
    if(options.poly):
      self.computeB = self.computeBPoly
    else:
      self.computeB = self.computeBLinear
      
    self.toleranceInner = options.maxToleranceInner
    
    if(options.useChannelWidth):
      self.computeW = self.computeWLinear
      
    Nx = options.Nx
    # Chebyshev operators
    self.cheb = Chebyshev.Chebyshev(Nx)
    
    self.dt = options.initDt
    
    # values for constants
    self.rho_w = 1000.0 # kg/m^3 water density
    self.sPerY = 365.25*24.*3600. # number of seconds per year
    self.g = 9.8 # m/s^2 gravity acceleration
    self.aBar = .3/self.sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m height scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    self.n = 3. # Glen's flow parameter
    self.WBar = 10000. # m width scaling factor
    self.tBar = self.xBar/self.uBar
    #print 'one year:', self.sPerY/self.tBar
    #exit()
    
    n = self.n
    rho_i = self.options.rho_i
    A = self.options.A
    Ab = self.options.Ab
    C = self.options.C
    lambda_0 = self.options.lambda_0
    m_0 = self.options.m_0
    NBar = rho_i*self.g*self.HBar
    tauLBar = A**(-1./n)*self.HBar*self.uBar**(1./n)/self.xBar**(1./n + 1)
    tauBBar = C*self.uBar**(1./n)
    tauDBar = rho_i*self.g*self.HBar**2/self.xBar
    self.Kappa = m_0/(lambda_0*Ab)*self.uBar/NBar**n
    self.epsilon = tauLBar/(2*tauDBar)
    self.gamma = tauBBar/tauDBar
    
    tauWBar = self.HBar/self.WBar*((n+2.)*self.uBar/(2.*A*self.WBar))**(1./n)
    self.omega = tauWBar/tauDBar
    self.delta = 1.0 - rho_i/self.rho_w

    self.outputFileIndex = 0
    self.time = 0

    if(options.inFile == "none" or options.inFile == "zero"):
      self.makeInitialCondition(options.inFile)
    else:
      if(not os.path.exists(options.inFile)):
        print "File not found:", options.inFile
        exit(1)
      self.readResults(options.inFile)

    self.oldH = self.H
    self.oldU = self.u
    self.oldXg = self.xg
    
    if(options.plot):
      self.prev_dH_dt = numpy.zeros(self.H.shape)
      self.prev_du_dt = numpy.zeros(self.u.shape)

    
    print "initial self:"
    pprint(self.__dict__)

  def makeInitialCondition(self, initOption):
    Nx = self.options.Nx
    a = self.options.a
    xc = self.options.xc
    sigma = self.cheb.x
    
    self.time = 0
    if(initOption == "none"):
      (HGuessSheet,xgGuess) = self.initH()
      if(xgGuess > self.options.xc):
        print "xgGuess=%f > xc=%f: cannot proceed!"%(xgGuess,self.options.xc)
        exit(1)
      HGuess = numpy.zeros((2*Nx))
      HGuess[0:Nx] = HGuessSheet
      uxg = a*xgGuess/HGuessSheet[-1]
      xShelf = xgGuess + (xc-xgGuess)*sigma
      const = (self.delta*a/(8.*self.epsilon))**self.n
      uShelf = (const*(xShelf**(self.n+1)-xgGuess**(self.n+1)) + uxg**(self.n+1))**(1./(self.n + 1))
      HGuess[Nx:2*Nx] = a*xShelf/uShelf
    
    if(initOption == "zero"):
      print "initializing with zero."
      HGuess = 1e-4*numpy.ones((2*Nx))
      xgGuess = 0.7
      
    
    self.updateHandXg(HGuess,xgGuess)
    self.u = numpy.zeros(self.H.shape)
    self.ux = numpy.zeros(self.H.shape)
    self.oldH = self.H
    self.oldU = self.u
    self.oldXg = self.xg
    self.dxg_dt = 0.
  
    innerConverged = False
    print "Computing initial u by iterating on viscosity:"
    for inner in range(self.options.maxInitSteps):
      prevU = self.u
      self.iterateOnViscosity()
      diffU = numpy.amax(numpy.abs(prevU-self.u))/numpy.amax(numpy.abs(self.u))
      if(self.options.plot):
        import matplotlib.pyplot as plt
        if(self.options.plotContinuous):
          plt.draw()
          plt.pause(0.0001)
        else:
          plt.show()
      print "iter: ", inner, "diffU:", diffU, "resStress: ", self.resStress
      if(diffU < self.options.initUTolerance):
        innerConverged = True
        break
  
    if not innerConverged:
      print "Warning: initial velocity did not converge after %i steps!"%self.options.maxInitSteps
      
    self.writeResults()
    
    cflSheet = numpy.amax(numpy.abs(self.u[0:Nx]/(self.xg*self.deltaSigma)))*self.dt
    cflShelf = numpy.amax(numpy.abs(self.u[Nx:2*Nx]/((xc-self.xg)*self.deltaSigma)))*self.dt
    cfl = numpy.maximum(cflSheet,cflShelf)

    # suggest a conservative estimate as the staring time step
    dt = self.dt*self.options.goalCFL/cfl
    print "Suggested dt:", dt, "or", dt*self.tBar/self.sPerY, "years"

  def initH(self):
  
    xg = scipy.optimize.brentq(bcRes, self.options.minXg, self.options.maxXg,
                               args=(self,))
    (b_xg,bx_xg) = self.computeB(xg)
    # initial conditions:
    H_xg = b_xg/(1-self.delta)
    
    y_xg = [H_xg]
    
    xInteg = self.cheb.x[::-1]*xg
    
    y = odeint(dy_dx, y_xg, xInteg, args=(self,))
  
    H = y[::-1,0]
    return (H,xg)
    
  def readResults(self, fileName):
    print "reading from: ", fileName
  
    filePointer = open(fileName,'rb')
    Nx = numpy.fromfile(filePointer, dtype=int, count=1)[0]
    if(Nx != self.options.Nx):
       print "Bad Nx in file."
       exit(1)
    
    self.xg = numpy.fromfile(filePointer, dtype=float, count=1)[0]
    self.time = numpy.fromfile(filePointer, dtype=float, count=1)[0]
    if self.writeToSeparateFiles and not self.initFromChebSteadyState:
       print "Initializing from restart file"
       self.outputFileIndex = numpy.fromfile(filePointer, dtype=int, count=1)[0]  
    else:
       print "Initializing from initFromChebSteadyState"
    self.H = numpy.fromfile(filePointer, dtype=float, count=2*Nx)
    self.updateHandXg(self.H,self.xg)
    self.u = numpy.fromfile(filePointer, dtype=float, count=2*Nx)
    self.ux = numpy.dot(self.Dx,self.u)
    self.dxg_dt = numpy.fromfile(filePointer, dtype=float, count=1)[0]
    filePointer.close()
  	
    solver.writeResults()
  
  def writeResults(self):
    if self.options.writeToSeparateFiles:
        parts = string.rsplit(self.options.outFile,'.',1)
        fileName = "%s.%04i.%s"%(parts[0],self.outputFileIndex,parts[1])
    else:
        fileName = "%s"%(self.options.outFile)
  
    print "writing to: ", fileName
    filePointer = open("%s/%s"%(self.options.folder,fileName),'wb')
    Nx = numpy.array(self.options.Nx,int)
    Nx.tofile(filePointer)
    xg = numpy.array(self.xg)
    xg.tofile(filePointer)
    time = numpy.array(self.time)
    time.tofile(filePointer)
    if self.options.writeToSeparateFiles:
       outputFileIndex = numpy.array(self.outputFileIndex,int)
       outputFileIndex.tofile(filePointer)
       print "outputFileIndex is %s"%(outputFileIndex)
    self.H.tofile(filePointer)
    self.u.tofile(filePointer)
    dxg_dt = numpy.array(self.dxg_dt)
    dxg_dt.tofile(filePointer)
    self.x.tofile(filePointer)
    self.sx.tofile(filePointer)
    self.ux.tofile(filePointer)
    p = numpy.array(self.options.p)
    p.tofile(filePointer)
    A = numpy.array(self.options.A)
    A.tofile(filePointer)
    C = numpy.array(self.options.C)
    C.tofile(filePointer)
    rho_i = numpy.array(self.options.rho_i)
    rho_i.tofile(filePointer)
    a = numpy.array(self.options.a)
    a.tofile(filePointer)
    meltRate = numpy.array(self.options.meltRate)
    meltRate.tofile(filePointer)
    linearSlope = numpy.array(self.options.linearSlope)
    linearSlope.tofile(filePointer)
    eps_s = numpy.array(self.options.eps_s)
    eps_s.tofile(filePointer)
    xc = numpy.array(self.options.xc)
    xc.tofile(filePointer)
    filePointer.close()
    filePointer = open("%s/%s"%(self.options.folder,self.options.filePointer),'w')
  #  filePointer.write("%s\n"%options.outFile)
    filePointer.write("%s\n"%(fileName))
    filePointer.close()
  
    self.outputFileIndex +=1

  # call each time H and xg are updated to compute b, Hf, Hx, bx and Np
  def updateHandXg(self, H, xg):
    Nx = self.options.Nx
    xc = self.options.xc
    p = self.options.p
    if(xg >= xc):
      print "Error: xg=%f >= xc=%f. Exiting..."%(xg,xc)
      exit(1)
    self.DxSheet = self.cheb.Dx/xg
    xSheet = xg*self.cheb.x
    self.DxShelf = self.cheb.Dx/(xc-xg)
    xShelf = xg + (xc-xg)*self.cheb.x
    
    self.deltaSigma = numpy.zeros(self.cheb.x.shape)
    self.deltaSigma[0:(Nx+1)/2] = self.cheb.x[1:(Nx+3)/2]-self.cheb.x[0:(Nx+1)/2]
    self.deltaSigma[(Nx+1)/2:] = self.cheb.x[(Nx+1)/2:]-self.cheb.x[(Nx-1)/2:-1]
    
    x = numpy.zeros((2*Nx))
    x[0:Nx] = xSheet
    x[Nx:2*Nx] = xShelf

    self.Dx = numpy.zeros((2*Nx,2*Nx))
    self.Dx[0:Nx,0:Nx] = self.DxSheet
    self.Dx[Nx:2*Nx,Nx:2*Nx] = self.DxShelf

    self.H = H
    self.xg = xg
    self.x = x
    
    self.Hx = numpy.dot(self.Dx,self.H)

    (self.b,bx) = self.computeB(x)
    self.Hf = self.maxS(self.b[0:Nx])/(1-self.delta)
    self.Hfx_xg = bx[Nx]/(1-self.delta)
    self.Np = self.H[0:Nx]*(self.maxS(1.0 - self.Hf/self.H[0:Nx]))**p

    self.s = numpy.zeros(x.shape)
    self.s[0:Nx] = self.H[0:Nx]-self.b[0:Nx]
    self.s[Nx:2*Nx] = self.delta*self.H[Nx:2*Nx]

    self.sx = numpy.zeros(x.shape)
    self.sx[0:Nx] = self.Hx[0:Nx]-bx[0:Nx]
    self.sx[Nx:2*Nx] = self.delta*self.Hx[Nx:2*Nx]
    
    if(self.options.useChannelWidth):
      (self.W,self.Wx) = self.computeW(x)
      
  def absS(self,x):
    return numpy.sqrt(x**2 + self.options.eps_s**2)
    
  def absSx(self,x):
    return x/self.absS(x)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)
    
    return 0.5*(self.absSx(x) + 1.0)
    
  def computeWLinear(self, x):
    W = self.options.W0 + x*self.options.Wx
    Wx = numpy.ones(x.shape)*self.options.Wx
    return (W,Wx)
  
  def computeBLinear(self,x):
    shoofx = 750000.
    slope = self.options.linearSlope/self.HBar
    point = 778.5/self.HBar
    intercept = -720/self.HBar
  
    bx = slope*self.xBar/shoofx
    b0 = slope*intercept/point
  
    eps_b = 1e-8
    absx = numpy.sqrt(x**2 + eps_b**2)
    b = b0 + absx*bx
    
    bx = bx*x/absx
   
    return (b,bx)
    
  def computeBPoly(self,x):
    schoofx = 750000.
    xs = x*self.xBar/schoofx
    b = -(729 - 2184.8*xs**2 + 1031.72*xs**4 - 151.72*xs**6)/self.HBar  
    bx = -self.xBar/schoofx*(- 2.0*2184.8*xs + 4.0*1031.72*xs**3 - 6.0*151.72*xs**5)/self.HBar  
    return (b,bx)
  
    
  def iterateOnViscosity(self):
    # make local copies of common variables for convenience
    Nx = self.options.Nx
    xc = self.options.xc
    x = self.x
    n = self.n
    
    # these variables should be computed by updateH()
    # before iteration begins
    H = self.H
    Hf = self.Hf
    Np = self.Np
    
    uk = self.u # previous iteration
    uxk = self.ux

    #tau_l = (H nu ux)_x
    #nu = 4*epsilon*|ux|^(1/n - 1)
    abs_ux = self.absS(uxk)
    tauLCoeff = 4.*self.epsilon*(abs_ux)**(1./n-1.)*H
    #if(not self.useLongi):
    #  tauLCoeff[0:Nx] = 0.

    #tau_b = -gamma*(N^n/(Kappa*|u| + N^n))^(1/n) u/|u|
    abs_u = self.absS(uk)
    tauBCoeff = numpy.zeros(x.shape)
    if(self.options.useSchoofBasal):
      tauBCoeff[0:Nx] = -self.gamma*abs_u[0:Nx]**(1/n-1)
    else:
      tauBCoeff[0:Nx] = -self.gamma*Np*(abs_u[0:Nx]/(self.Kappa*abs_u[0:Nx] + Np**n))**(1./n)/abs_u[0:Nx]

    #tau_w = -omega*H*W^(-1-1/n)*|u|^(1/n-1) u
    if(self.options.useChannelWidth):
      tauWCoeff = -self.omega*H*self.W**(-1.-1./n)*abs_u**(1/n-1)
    else:
      tauWCoeff = numpy.zeros(x.shape)

    #tau_d = -H sx
    tauD = -H*self.sx

    tauLk = numpy.dot(self.Dx,tauLCoeff*uxk)
    tauBk = tauBCoeff*uk
    tauWk = tauWCoeff*uk
    
    residual = tauLk + tauBk + tauD + tauWk
    residual = tauLk + tauBk + tauD + tauWk
    
    M = numpy.dot(self.Dx,numpy.dot(numpy.diag(tauLCoeff),self.Dx)) \
      + numpy.diag(tauBCoeff + tauWCoeff)
    rhs = -tauD
    
    # velocity is zero at the ice divide
    M[0,:] = 0.
    M[0,0] = 1.
    rhs[0] =  0.
    
    # H nu ux = 0.5*delta*H^2 at x = xc
    tauLCoeff_xg = 4.*self.epsilon*(abs_ux[-1])**(1./n-1.)*H[-1]
    M[-1,:] = self.Dx[-1,:]
    rhs[-1] = 0.5*self.delta*H[-1]**2/tauLCoeff_xg
    
    # have to replace rows to get continuity in u and ux at xg
    M[Nx-1,:] = 0.
    M[Nx-1,Nx-1] = 1.
    M[Nx-1,Nx] = -1.
    rhs[Nx-1] = 0.
    M[Nx,0:Nx] = self.DxSheet[Nx-1,:]
    M[Nx,Nx:2*Nx] = -self.DxShelf[0,:]
    rhs[Nx] = 0.
        
    ukp1 = linalg.solve(M,rhs)
    uxkp1 = numpy.dot(self.Dx,ukp1)
    
    #self.noiseFloor = numpy.max(numpy.abs(self_res[1:-1]))

    if(self.options.plot):
      import matplotlib.pyplot as plt
      residual_check = numpy.dot(M,uk) - rhs
      self_res = numpy.dot(M,ukp1) - rhs
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,1)
      ax.cla()
      oldXSheet = self.cheb.x*self.oldXg
      oldXShelf = self.oldXg + (xc-self.oldXg)*self.cheb.x
      (oldB,oldBx) = self.computeB(oldXSheet)
      plt.plot(x, -self.b, 'k', x[0:Nx], Hf-self.b[0:Nx], 'g',
               x[Nx:2*Nx],-(1.-self.delta)*self.H[Nx:2*Nx], 'r', 
               oldXSheet,self.oldH[0:Nx]-oldB,'b',
               oldXShelf,self.delta*self.oldH[Nx:2*Nx],'b',
               x,self.s, 'r')
      ax = plt.subplot(2,3,2)
      ax.cla()
      plt.plot(x[1:-1],self_res[1:-1],'b', x[1:-1], residual[1:-1], 'r',
               x[1:-1], residual_check[1:-1], 'k--')
      ax = plt.subplot(2,3,3)
      ax.cla()
      plt.plot(x,tauLk, 'b', x, tauBk, 'r', x, tauWk, 'm', x, tauD, 'g', x, residual, 'k')
      ax = plt.subplot(2,3,4)
      ax.cla()
      plt.plot(x,uk, 'b', x, ukp1, 'r', x[0:Nx], Np**n/self.Kappa,'g')
      plt.ylim(numpy.amin(ukp1),numpy.amax(ukp1))
      
      
    self.resStress = numpy.max(numpy.abs(residual[1:-1]))
    
    self.u = ukp1
    self.ux = uxkp1
    
  def newtonStepXg(self, newH):
    Nx = self.options.Nx
    # Newton's method to find the new xg
    # xg^(k+1) = xg^k - f(xg^k)/f'(xg^k)
    # f = Hf-H
    f = self.Hf - newH[0:Nx]
    fx_xg = numpy.dot(self.Dx[Nx-1,0:Nx],f)
    newXg = self.xg - f[-1]/fx_xg
    
    return newXg

  def compute_dxg_dt(self):    
    Nx = self.options.Nx
    a = self.a
    u = self.u[Nx]
    H = self.H[Nx]
    ux = self.ux[Nx]
    Hx = self.Hx[Nx]
    Wx = self.Wx[Nx]
    W = self.W[Nx]
    Hfx_xg = self.Hfx_xg
    
    # dH_dt = dHf_dt
    # a - (u - dxg_dt) Hx - ux H - (u - dxg_dt) H Wx/W = Hfx*dxg_dt
    # dxg_dt = -(a - u Hx - ux H - Hu Wx/W)/(Hx + H Wx/W - Hfx)
    numerator = a - u*Hx - ux*H
    denominator = Hx - Hfx_xg
    if(self.options.useChannelWidth):
      numerator -= H*u*Wx/W
      denominator += H*Wx/W
    
    dxg_dt = -numerator/denominator
    return dxg_dt
    
    
  def computeExplicit_dH_dt(self,dxg_dt): 
    Nx = self.options.Nx
    a = self.options.a
    u = self.u
    H = self.H
    ux = self.ux
    Hx = self.Hx
    sigma = self.cheb.x
    
    # dH/dt + (u - sigma dxg_dt) Hx + ux H + (u - sigma dxg_dt) H Wx/W = a (in sheet)
   
    movingGridTerm = numpy.zeros(self.x.shape)
    movingGridTerm[0:Nx] = -sigma*dxg_dt
    movingGridTerm[Nx:2*Nx] = -(1.-sigma)*dxg_dt
   
    meltRate = numpy.zeros(self.x.shape)
    meltRate[Nx:2*Nx] = self.options.meltRate
    
    dH_dt = a - meltRate -(u + movingGridTerm)*Hx - ux*H
    if(self.options.useChannelWidth):
      dH_dt -= (u + movingGridTerm)*H*self.Wx/self.W
    
    return dH_dt
    
  def iterateImplicitContinuity(self):
    # u Hx + ux H = a
    sigma = self.cheb.x
    Nx = self.options.Nx 
    
    dxg_dt = self.dxg_dt
        
    movingGridTerm = numpy.zeros(self.x.shape)
    movingGridTerm[0:Nx] = -sigma*dxg_dt
    movingGridTerm[Nx:2*Nx] = -(1.-sigma)*dxg_dt

    meltRate = numpy.zeros(self.x.shape)
    meltRate[Nx:2*Nx] = self.options.meltRate
    
    theta = self.options.timeCentering   
   
    diag1 = 1./self.dt + theta*self.ux
    if(self.options.useChannelWidth):
      diag1 += theta*(self.u + movingGridTerm)*self.Wx/self.W
    

    # dH/dt + (u - sigma dxg_dt) Hx + ux H + (u - sigma dxg_dt) H Wx/W = a (in sheet)
    M =  numpy.diag(diag1) + numpy.dot(numpy.diag(theta*(self.u + movingGridTerm)),self.Dx)
    rhs = theta*(self.options.a - meltRate) + self.oldH/self.dt + (1. - theta)*self.old_dH_dt
    
    # have to replace rows to get continuity in H and Hx at xg
    M[Nx-1,:] = 0.
    M[Nx-1,Nx-1] = 1.
    M[Nx-1,Nx] = -1.
    rhs[Nx-1] = 0.
    M[Nx,0:Nx] = self.DxSheet[Nx-1,:]
    M[Nx,Nx:2*Nx] = -self.DxShelf[0,:]
    rhs[Nx] = 0.
    
    # no boundary conditions should be necessary -- ux H = a at the ice divide
    newH = linalg.solve(M,rhs)
    
    if(self.options.plot):
      import matplotlib.pyplot as plt
      x = self.x
      newHx = numpy.dot(self.Dx,newH)
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,6)
      ax.cla()
      plt.plot(x,(newH-self.oldH)/self.dt, 'b', x,movingGridTerm*newHx, 'r', 
               x, self.u*newHx, 'k', x, self.ux*newH, 'g')
    return newH

  def takeTimeStep(self):
    Nx = self.options.Nx 
    xc = self.options.xc
    sigma = self.cheb.x
    
    self.oldH = self.H
    self.oldU = self.u
    self.oldXg = self.xg
    print self.dxg_dt
    self.old_dH_dt = self.computeExplicit_dH_dt(self.dxg_dt)
    innerConverged = False
    
    # initial guess is that du/dt will be the same over this step
    newH = self.iterateImplicitContinuity()
    newXg = self.newtonStepXg(newH)
    self.updateHandXg(newH,newXg)
  
    for inner in range(self.options.maxInnerSteps):
      prevU = self.u
      prevH = self.H
      self.iterateOnViscosity()
      newH = self.iterateImplicitContinuity()
      newXg = self.newtonStepXg(newH)
  
      diffH = numpy.amax(numpy.abs(newH-prevH))/numpy.amax(newH)
      diffXg = numpy.abs(newXg-self.xg)/newXg
      self.updateHandXg(newH,newXg)
  
      diffU = numpy.amax(numpy.abs(prevU-self.u))/numpy.amax(numpy.abs(self.u))
  
      self.dxg_dt = (newXg- self.oldXg)/self.dt
      
      uEff = self.u[0:Nx] - sigma*self.dxg_dt
      cflSheet = numpy.amax(numpy.abs(uEff/(newXg*self.deltaSigma)))*self.dt
      uEff = self.u[Nx:2*Nx] - (1.-sigma)*self.dxg_dt
      cflShelf = numpy.amax(numpy.abs(uEff/((xc-newXg)*self.deltaSigma)))*self.dt
      cfl = numpy.maximum(cflSheet,cflShelf)
      
      print cflSheet, cflShelf
          
      if numpy.isnan(cfl):
        print "blew up!"
        exit(1)
        
      print "time: ", self.time+self.dt, "iter: ", inner, "cfl:", cfl, "diffs: ", diffH, diffXg, diffU, self.resStress
      
      if(self.options.plot):
        import matplotlib.pyplot as plt
        dH_dt = (newH-self.oldH)/self.dt
        du_dt = (self.u-self.oldU)/self.dt
        
        xSheet = self.xg*sigma
        xShelf = self.xg + (xc-self.xg)*sigma
        x = numpy.zeros((2*Nx))
        x[0:Nx] = xSheet
        x[Nx:2*Nx] = xShelf
  
        ax = plt.subplot(2,3,5)
        ax.cla()
        plt.plot(x,self.prev_dH_dt, 'g')
        plt.plot(x,self.prev_du_dt, 'm')
        plt.plot(x,dH_dt, 'b')
        plt.plot(x,du_dt, 'r')
  
        dH_dt = numpy.amax(numpy.abs(dH_dt))
        ax = plt.subplot(2,3,1)
        plt.title('iter: %02i dt=%.4g'%(inner, self.dt))
        ax = plt.subplot(2,3,2)
        plt.title('CFL=%.2f'%cfl)
        ax = plt.subplot(2,3,3)
        plt.title('diffU: %.4g tol.: %.4g'%(diffU, self.toleranceInner))
        ax = plt.subplot(2,3,4)
        plt.title('xg=%.4f'%(newXg))
        ax = plt.subplot(2,3,5)
        plt.title('dxg/dt=%.4f'%(self.dxg_dt))
        ax = plt.subplot(2,3,6)
        plt.title('|dH/dt|_max=%.4f'%(dH_dt))
        if(self.options.plotContinuous):
          plt.draw()
          plt.pause(0.0001)
        else:
          plt.show()
  
      if(diffH < self.toleranceInner and diffXg < self.toleranceInner and diffU < self.toleranceInner):
        innerConverged = True
        break
      
    if not innerConverged:
      print "Error: inner loop did not converge after %i steps!"%self.options.maxInnerSteps
      print "Try reducing the goal CFL number."
      exit(1)
    
    dH_dt = (newH-self.oldH)/self.dt
    if(self.options.plot):
      self.prev_dH_dt = dH_dt
      self.prev_du_dt = (self.u-self.oldU)/self.dt
    #self.dxg_dt = (newXg-self.oldXg)/self.dt
    diffH = numpy.amax(numpy.abs(dH_dt))
    diffXg = numpy.abs(self.dxg_dt)
    self.time += self.dt
    if(diffH < self.options.toleranceH and diffXg < self.options.toleranceXg):
      return (True,True)
    # make the inner tolerance stricter as H and xg converge
    # changes in the inner loop should be at most the size of changes in the outer loop,
    # once we have reached the goal CFL
    maxChange = max(diffH/numpy.amax(newH),diffXg/newXg)
    self.toleranceInner = min(self.options.maxToleranceInner,self.dt*self.options.goalCFL/cfl*maxChange)
  
    if not self.options.fixedTimeStep:
        scale = max(0.5,min(2.0,self.goalCFL/cfl))
        self.dt *= scale
  
    print "time: ", self.time, "|dH_dt|_max: ", diffH, "dxg_dt:", self.dxg_dt, "dt: ", self.dt, "inner tol.:", self.toleranceInner
    self.updateHandXg(newH,newXg)
   
    if self.time>=self.options.finalTime/self.tBar*self.sPerY:
      print "The run exceeded the transient final time."
      return (True,False)
   
    return (False,False)
      

parser = OptionParser()

parser.add_option("--xc", type="float", default=2.112, dest="xc")
parser.add_option("--useChannelWidth", action="store_true", dest="useChannelWidth")
parser.add_option("--W0", type="float", default=1000, dest="W0")
parser.add_option("--Wx", type="float", default=0.0, dest="Wx")

parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=1e-25, dest="A")
parser.add_option("--Ab", type="float", default=3.1688e-24, dest="Ab")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--rho_i", type="float", default=900.0, dest="rho_i")
parser.add_option("--a", type="float", default=1.0, dest="a")
parser.add_option("--meltRate", type="float", default=0.0, dest="meltRate")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007
parser.add_option("--lambda_0", type="float", default=2, dest="lambda_0")
parser.add_option("--m_0", type="float", default=0.5, dest="m_0")
parser.add_option("--poly", action="store_true", dest="poly")

parser.add_option("--inFile", type="string", default="none", dest="inFile")
parser.add_option("--folder", type="string", default="results", dest="folder")
parser.add_option("--outFile", type="string", default="results.pyda", dest="outFile")
parser.add_option("--filePointer", type="string", default="default.pointer", dest="filePointer")

parser.add_option("--Nx", type="int", default=513, dest="Nx")
parser.add_option("--maxSteps", type="int", default=1000, dest="maxSteps")
parser.add_option("--maxInnerSteps", type="int", default=200, dest="maxInnerSteps")
parser.add_option("--maxInitSteps", type="int", default=200, dest="maxInitSteps")
parser.add_option("--stepsPerWrite", type="int", default=10, dest="stepsPerWrite")
parser.add_option("--timeCentering", type="float", default=1.0, dest="timeCentering")
parser.add_option("--writeToSeparateFiles", action="store_true", default=False, dest="writeToSeparateFiles")
parser.add_option("--initFromChebSteadyState", action="store_true", default=False, dest="initFromChebSteadyState")
parser.add_option("--finalTime", type="float", default=500, dest="finalTime")
parser.add_option("--fixedTimeStep", action="store_true", default=False, dest="fixedTimeStep")
parser.add_option("--useSchoofBasal", action="store_true", default=False, dest="useSchoofBasal")

parser.add_option("--goalCFL", type="float", default=2.0, dest="goalCFL")
parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")
parser.add_option("--toleranceH", type="float", default=0.1, dest="toleranceH")  # 0.03 m/yr
parser.add_option("--toleranceXg", type="float", default=0.1, dest="toleranceXg")  # 30 m/yr
parser.add_option("--maxToleranceInner", type="float", default=1e-3, dest="maxToleranceInner")
parser.add_option("--initUTolerance", type="float", default=1e-3, dest="initUTolerance")
parser.add_option("--initDt", type="float", default=1.5e-4, dest="initDt")
parser.add_option("--minXg", type="float", default=0.75, dest="minXg")
parser.add_option("--maxXg", type="float", default=2.0, dest="maxXg")

parser.add_option("--plot", action="store_true", dest="plot")
parser.add_option("--plotContinuous", action="store_true", dest="plotContinuous")

optionsStruct, args = parser.parse_args()

print "options:"
pprint(optionsStruct.__dict__)

solver = SheetShelfSolver(optionsStruct)

for outer in range(solver.options.maxSteps):
  (finished,converged) = solver.takeTimeStep()
  if(finished):
    break
  if(numpy.mod(outer,solver.options.stepsPerWrite) == solver.options.stepsPerWrite-1):
    solver.writeResults()
    
if(converged):
  print "The run converged."
elif(not finished):
  print "The run did not converge after %i steps."%solver.options.maxSteps

solver.writeResults()

if(not finished):
  exit(1)
