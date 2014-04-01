from pprint import pprint
import numpy
import scipy.sparse
import scipy.sparse.linalg
import numpy.linalg
import os.path
#import matplotlib.pyplot as plt

class Solver:
  def __init__(self, options):
    if(options.poly):
      self.computeB = self.computeBPoly
    else:
      self.computeB = self.computeBLinear
    self.Nx = options.Nx
    self.deltaX = options.deltaX
    self.dt = options.dt
    self.xc = options.xc
    self.xu = numpy.linspace(0,self.xc,self.Nx)
    self.xH = self.xu-0.5*self.deltaX
    self.rho_i = options.rho_i 
    self.rho_w = options.rho_w 
    self.delta = 1.0 - self.rho_i/self.rho_w
    self.C = options.C
    self.W = options.W
    self.a = options.a
    self.p = options.p
    self.A = options.A
    self.g = options.g # m/s^2 gravity acceleration
    self.n = options.n # Glen's flow parameter
    self.linearSlope = options.linearSlope

    self.Ab = options.Ab
    self.lambda_0 = options.lambda_0 # wavelength of bedrock bump (m)
    self.m_0  = options.m_0 # maximum bed obstacle slope (no unit)

    self.plotContinuous = options.plotContinuous
    self.plot = options.plotContinuous
    self.disableLongi = options.disableLongi
    self.useChannel = options.useChannel
    self.useSchoofBasal = options.useSchoofBasal
    self.useGLP = options.useGLP
    self.maxPicardIter = options.maxInnerSteps
    
    self.sPerY = 365.25*24.*3600. # number of seconds per year
    self.aBar = .3/self.sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.WBar = 10000. # m channel width scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor    
    self.tBar = self.xBar/self.uBar    
    
    self.eps_s = options.eps_s
    self.toleranceInner = options.toleranceInner
    self.tol_out_dH_dt = options.toleranceH/self.sPerY/self.aBar # tolerance equivaent to 0.1 m/yr
    self.tol_out_dxg_dt = options.toleranceXg/self.sPerY/self.uBar # tolerance equivaent to 1 mm/yr (MISMIP)

    self.NBar = self.rho_i*self.g*self.HBar
    self.taudBar = self.NBar*self.HBar/self.xBar
    self.taubBar = self.C*self.uBar**(1./self.n)
    self.Kappa = (self.m_0*self.uBar)/(self.lambda_0*self.Ab)/self.NBar**self.n
    self.gamma = self.taubBar/self.taudBar

    self.taulBar = (self.A**(-1./self.n))*(self.uBar/self.xBar)**(1./self.n)*self.HBar/self.xBar
    self.epsilon = self.taulBar/(2*self.taudBar)

    self.tauWbar = self.HBar/self.WBar*(5*self.uBar/(2*self.A*self.WBar))**(1./self.n)
    self.omega = self.tauWbar/self.taudBar
    
    self.outFile = options.outFile
    self.folder = options.folder
    self.filePointer = options.filePointer
    #print 'gamma, Kappa, epsilon, omega:', self.gamma, self.Kappa, self.epsilon, self.omega


    diagonals = []
    diag0 = -1./self.deltaX*numpy.ones((self.Nx),float)
    diag0[-1] = 0
    diagonals.append(diag0)
    diag1 = 1./self.deltaX*numpy.ones((self.Nx-1),float)
    diagonals.append(diag1)
    self.DxH = scipy.sparse.diags(diagonals, [0,1], shape=(self.Nx,self.Nx), format='csr')
    
    diagonals = []
    diag0 = 1./self.deltaX*numpy.ones((self.Nx),float)
    diag0[0] = 0
    diagonals.append(diag0)
    diag1 = -1./self.deltaX*numpy.ones((self.Nx-1),float)
    diagonals.append(diag1)
    self.Dxu = scipy.sparse.diags(diagonals, [0,-1], shape=(self.Nx,self.Nx), format='csr')
 
    diagonals = []
    diag0 = 0.5*numpy.ones((self.Nx),float)
    diag0[-1] = 0
    diagonals.append(diag0)
    diag1 = 0.5*numpy.ones((self.Nx-1),float)
    diagonals.append(diag1)
    self.AvgH = scipy.sparse.diags(diagonals, [0,1], shape=(self.Nx,self.Nx), format='csr')
    
    diagonals = []
    diag0 = 0.5*numpy.ones((self.Nx),float)
    diag0[0] = 0
    diagonals.append(diag0)
    diag1 = 0.5*numpy.ones((self.Nx-1),float)
    diagonals.append(diag1)
    self.Avgu = scipy.sparse.diags(diagonals, [0,-1], shape=(self.Nx,self.Nx), format='csr')

    (self.b,self.bx) = self.computeB(self.xH)
    #print 'b:', self.b
    self.Hf = self.maxS(self.b)/(1-self.delta)
    self.bx = self.DxH.dot(self.b)
    
    print "options:"
    pprint(options.__dict__)
    
    print "self:"
    pprint(self.__dict__)    
    
    if(self.plot):
      import matplotlib.pyplot as plt

    if(self.plot and self.plotContinuous):
      plt.ion()
      
    if(options.inFile == "none" or options.inFile == "zero"):
      self.makeInitialCondition(options.inFile, options.xgInit, options.maxInitSteps,options.initUTolerance)
    else:
      if(not os.path.exists(options.inFile)):
        print "File not found:", options.inFile
        exit(1)
      self.readResults(options.inFile)

  def makeInitialCondition(self, inFile, xgInit, maxInitSteps, initUTolerance):
    self.time = 0
    if(inFile == "none"):
      (HGuess,uGuess) = self.computeHGuess(xgInit)
    else:
      HGuess = 1e-4*numpy.ones(self.xH)
      uGuess = 1e-4*numpy.ones(self.xu)
  
    self.H = HGuess
    self.u = uGuess
    self.updateXgAndFlotationMasks()
    self.HPrev = self.H
    self.uPrev = self.u
    self.xgPrev = self.xg

   # Check (reinitializing) uGuess using Picard iteration 
    innerConverged = False
    print "Computing initial u by iterating on viscosity:"
    for inner in range(maxInitSteps):
      prevU = self.u
      self.iterateOnViscosity()
      diffU = numpy.amax(numpy.abs(prevU-self.u))/numpy.amax(numpy.abs(self.u))
      self.plotIter(0,inner)
      print "iter: ", inner, "diffU:", diffU, "resStress: ", self.resStress
      if(diffU < initUTolerance):
        innerConverged = True
        break

    self.innerConverged = innerConverged

    if not innerConverged:
      print "Warning: initial velocity did not converge after %i steps!"%maxInitSteps
    
    self.writeResults()
    
  def computeHGuess(self,xg):
    xH = self.xH
    xu = self.xu
    n = self.n
    a = self.a
    
    glIndex = numpy.argmin(numpy.abs(xH-xg))
    xg = xH[glIndex]
    #shelfMask = xH > xg
    
    (bxg,bxxg) = self.computeB(xg)
    Hxg = bxg/(1.0-self.delta)
    if(Hxg <= 0.):
      raise ValueError("Hfself(xg) <= 0. Cannot produce HGuess")
        
    (b,bx) = self.computeB(xH)
    
    uxg = a*xg/Hxg
    c1 = (self.delta*a/(8*self.epsilon))**n
    
    operand = numpy.maximum(c1*(xH**(n+1) - xg**(n+1)),0.0) + uxg**(n+1)
    uShelf = (operand)**(1/(n+1))
    HShelf = a*xH/uShelf
    
    H = HShelf.copy()
    u = uShelf.copy()
    
    # assume balance between taub and taud in the sheet (SIA) and steady state (u=a*x/H),
    # and iteratively solve the ODE for H
    Hip1 = H[glIndex+1]
    for xIndex in range(glIndex,-1,-1):
      Hi = Hip1
      deltaB = (b[glIndex+1]-b[glIndex])
      for iterIndex in range(100):
        Hmid = 0.5*(Hi+Hip1)
        umid = a*xu[xIndex]/Hmid
        taub = -self.gamma*umid**(1./n)
        #taud = -0.5*(Hi+Hip1)*((Hip1-Hi)/deltaX - deltaB/deltaX)
        # Hi**2 + Hi*deltaB + (-Hip1**2 +Hip1*deltaB + 2*deltaX*taub) = 0
        HPrev = Hi
        Hi = 0.5*(-deltaB + numpy.sqrt(deltaB**2-4*(-Hip1**2 +Hip1*deltaB + 2*self.deltaX*taub)))
        deltaH = numpy.abs(Hi-HPrev)
        if(deltaH < self.toleranceInner):
          break
      #print "deltaH:", Hi-HPrev, Hi, Hip1
      Hip1 = Hi
      H[xIndex] = Hi
      u[xIndex] = umid
    
  
    return (H,u)
    
  def readResults(self, fileName):
    print "reading from: ", fileName
    filePointer = open(fileName,'rb')
    Nx = numpy.fromfile(filePointer, dtype=int, count=1)[0]
    if(Nx != self.Nx):
      print "Bad Nx in file."
      exit(1)
  
    self.xg = numpy.fromfile(filePointer, dtype=float, count=1)[0]
    self.time = numpy.fromfile(filePointer, dtype=float, count=1)[0]
    self.H = numpy.fromfile(filePointer, dtype=float, count=Nx)
    self.updateXgAndFlotationMasks()
    self.u = numpy.fromfile(filePointer, dtype=float, count=Nx)
    filePointer.close()

  def writeResults(self):
    fileName = "%s/%s"%(self.folder,self.outFile)
    print "writing to: ", fileName
    filePointer = open(fileName,'wb')
    Nx = numpy.array(self.Nx,int)
    Nx.tofile(filePointer)
    xg = numpy.array(self.xg)
    xg.tofile(filePointer)
    time = numpy.array(self.time)
    time.tofile(filePointer)
    self.H.tofile(filePointer)
    self.u.tofile(filePointer)
    xH = self.xH
    xH.tofile(filePointer)
    xu = self.xu
    xu.tofile(filePointer)
  #  self.Hx.tofile(filePointer)
  #  self.ux.tofile(filePointer)
    p = numpy.array(self.p)
    p.tofile(filePointer)
    A = numpy.array(self.A)
    A.tofile(filePointer)
    C = numpy.array(self.C)
    C.tofile(filePointer)
    rho_i = numpy.array(self.rho_i)
    rho_i.tofile(filePointer)
    a = numpy.array(self.a)
    a.tofile(filePointer)
    linearSlope = numpy.array(self.linearSlope)
    linearSlope.tofile(filePointer)
    eps_s = numpy.array(self.eps_s)
    eps_s.tofile(filePointer)
    filePointer.close()
    filePointer = open(self.filePointer,'w')
    filePointer.write("%s\n"%self.outFile)
    filePointer.close()

  def computeBLinear(self,x):
      shoofx = 750000.
      slope = self.linearSlope/self.HBar
      pente = 778.5/self.HBar
      intercept = -720/self.HBar

      bx = self.linearSlope*self.xBar/shoofx/self.HBar
      b0 = slope*intercept/pente
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


  def absS(self,x):
    return numpy.sqrt(x**2 + self.eps_s**2)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)


#  def avgH(self,H):
#    H_u = numpy.zeros(H.shape,float)
#    H_u[0:-1] = 0.5*(H[1:]+H[0:-1])
#    return H_u
#
#  def avgu(self,u):
#    u_H = numpy.zeros(u.shape,float)
#    u_H[1:] = 0.5*(u[1:]+u[0:-1])
#    return u_H

  def plotIter(self,outer,inner):
    if(not self.plot):
      return
      
    import matplotlib.pyplot as plt
    if(self.plotContinuous):
      plt.ion()
        
    plt.figure(1, figsize=(16, 12))
    dH_dt = (self.H-self.HPrev)/self.dt
    dxg_dt = (self.xg-self.xgPrev)/self.dt
    du_dt = (self.u-self.uPrev)/self.dt
    
    ax = plt.subplot(2,3,1)
    ax.cla()
    s = self.H-self.b
    s[self.floatingMaskH]=self.delta*self.H[self.floatingMaskH]
    sPrev = self.HPrev-self.b
    mask = self.HPrev < self.Hf
    sPrev[mask] = self.delta*self.HPrev[mask]
    plt.plot(self.xH,s, 'r', self.xH, -self.b, 'k', 
             self.xH, self.Hf-self.b, 'g',self.xH, sPrev,'b')
    plt.legend(('s$^{t+1}$','b','H_f-b','s$^t$'))
             
    ax = plt.subplot(2,3,2)
    ax.cla()
    plt.plot(self.xu[1:-1], self.resTau_solver[1:-1],'b', 
             self.xu[1:-1], self.resTau_k[1:-1], 'r')
    plt.legend(('resTau$_{solver}$','resTau$^k$'))

    ax = plt.subplot(2,3,3)
    ax.cla()
    plt.plot(self.xu, self.longi, 'b', self.xu, self.basal, 'r', 
             self.xu, self.driving, 'g', self.xu, self.resTau_k, 'k')
    plt.legend(('longi','basal','driving','resTau$^k$'))

    ax = plt.subplot(2,3,4)
    ax.cla()
    plt.plot(self.xu,self.u, 'r', self.xu, self.uPrev, 'b')
    plt.ylim(numpy.amin(self.u),numpy.amax(self.u))
    plt.legend(('u$^{t+1}$','u$^t$'))

    ax = plt.subplot(2,3,5)
    ax.cla()
    plt.plot(self.xH,dH_dt, 'b')
    plt.plot(self.xu,du_dt, 'r')
    plt.legend(('dH/dt','du/dt'))
    
    ax = plt.subplot(2,3,6)
    ax.cla()
    plt.plot(self.xH, dH_dt, 'b', 
             self.xH, self.Dxu.dot(self.H*self.u), 'k', 
             self.xH, self.a*numpy.ones(self.xH.shape), 'r')
    plt.legend(('dH/dt','(uH)$_x$','a'))
   

    cfl = numpy.amax(numpy.abs(self.u/self.deltaX))*self.dt
    plt.subplot(2,3,1)
    plt.title('outer: %i inner: %i time: %g'%(outer, inner,self.time))
    plt.subplot(2,3,2)
    plt.title('CFL=%.2f'%cfl)
    plt.subplot(2,3,3)
    plt.title('resStress: %.4g tol.: %.4g'%(self.resStress, self.toleranceInner))
    plt.subplot(2,3,4)
    plt.title('lambda_g=%.4g xg=%.6g'%(self.lambda_g, self.xg))
    plt.subplot(2,3,5)
    plt.title('dxg/dt=%.4g tol.: %.4g'%(dxg_dt, self.tol_out_dxg_dt))
    plt.subplot(2,3,6)
    plt.title('|dH/dt|=%.4g tol.: %.4g'%(numpy.amax(numpy.abs(dH_dt)),self.tol_out_dH_dt))
    if(self.plotContinuous):
      plt.draw()
    else:
      plt.show()

        
        
  def computeSx(self):
    s = self.H-self.b
    s_floating = self.delta*self.H
    s[self.floatingMaskH] = s_floating[self.floatingMaskH]
    
    sx = self.DxH.dot(s)
    sx_floating = self.DxH.dot(s_floating)
    sx[self.floatingMaskU] = sx_floating[self.floatingMaskU]
    
    if(self.useGLP):
        sx[self.glIndices] = self.lambda_g*sx[self.glIndices] \
        + (1. - self.lambda_g)*sx_floating[self.glIndices]
    
    self.s = s
    self.sx = sx
    
    
  def updateXgAndFlotationMasks(self):
    fPattyn = self.Hf/self.H
    
    self.groundedMaskH = fPattyn < 1.
    self.floatingMaskH = fPattyn >= 1. 
    self.groundedMaskU = numpy.zeros(self.groundedMaskH.shape,bool)
    self.groundedMaskU[0:-1] = numpy.logical_and(self.groundedMaskH[1:],self.groundedMaskH[0:-1])
    
    self.floatingMaskU = numpy.ones(self.floatingMaskH.shape,bool)
    self.floatingMaskU[0:-1] = numpy.logical_and(self.floatingMaskH[1:],self.floatingMaskH[0:-1])
    
    self.glIndices = numpy.nonzero(numpy.logical_xor(self.groundedMaskH[1:],self.groundedMaskH[0:-1]))[0]
    if(len(self.glIndices) == 1):
      self.glIndices = self.glIndices[0]
    #print 'gl at index (or indices):', glIndices
        
    self.lambda_g = (1. - fPattyn[self.glIndices])/(fPattyn[self.glIndices+1] - fPattyn[self.glIndices])
    self.xg = self.lambda_g*self.xH[self.glIndices+1]+(1-self.lambda_g)*self.xH[self.glIndices]


  def iterateOnViscosity(self):
    p = self.p
    n = self.n # Glen's flow parameter
    #xH = self.xH
    #xu = self.xu
    Hk = self.H
    uk = self.u
    
    
    uxk = self.Dxu.dot(uk)
    abs_ux = self.absS(uxk)
    nuk = self.epsilon*(abs_ux)**(1./n-1.)
    
  
    # setting up tau_b  
     
    abs_u = self.absS(uk)
    if(self.useSchoofBasal):
      basalCoeff = -self.gamma*abs_u**(1./n-1)
    else:
      Nk = Hk*(self.maxS(1.0 - self.Hf/Hk))**p
      self.N = Nk
      Nk_u = self.AvgH.dot(Nk)
      Nk_u[self.floatingMaskU] = 0.
      #Nk_u[glIndices] = numpy.maximum(Nk[glIndices],Nk[glIndices+1])
      self.Nu = Nk_u
      gammaEff = self.gamma*(Nk_u**n/(self.Kappa*abs_u + Nk_u**n))**(1./n)
      self.basalu = -self.gamma*abs_u**(1./n-1)*uk
      self.basalN = -self.gamma/(self.Kappa)**(1/n)*Nk_u*abs_u**(-1)*uk
      basalCoeff = -gammaEff*abs_u**(1./n-1)
                 
    #print self.gamma
    if(self.useGLP):
        basalCoeff[self.glIndices] *= self.lambda_g #*2/(1+p)
   
    basalCoeff[self.floatingMaskU] = 0.
  
  
    # setting up Tau_l        
    
    longiCoeff = 4*Hk*nuk
    if(self.disableLongi):
      longiCoeff[self.groundedMaskU] = 0.
      

    # setting up Tau_w

    if(not self.useChannel):
        lateralEff = 0.
    else:
        lateralEff = self.omega*Hk*self.W**(1./n+1)
      
    lateralCoeff = -lateralEff*abs_u**(1./n-1)

  
    # setting up Tau_d
  
    self.computeSx()      
    Hk_u = self.AvgH.dot(Hk)
    driving_k = -Hk_u*self.sx
    self.driving = driving_k
    
    # Boundary condition at the calving front
    calving_ux = (self.delta/(8*self.epsilon)*Hk[-1])**n    
    
    # Residual calculation
    
    longi_k = self.DxH.dot(longiCoeff*uxk)

    self.longi = longi_k
    basal_k = basalCoeff*uk
    self.basal = basal_k
    lateral_k = lateralCoeff*uk
    self.lateral = lateral_k

    self.resTau_k = longi_k + basal_k + driving_k + lateral_k
    self.resTau_k[0] = uk[0]
    self.resTau_k[-1] = uxk[-1] - calving_ux
    
    # Setting up the matrix system to solve

    Mlon = scipy.sparse.diags([longiCoeff],[0],shape=(self.Nx,self.Nx), format='csr')
    Mlon = self.DxH.dot(Mlon.dot(self.Dxu))
    Mbas = scipy.sparse.diags([basalCoeff],[0],shape=(self.Nx,self.Nx), format='csr')
    Mlat = scipy.sparse.diags([lateralCoeff],[0],shape=(self.Nx,self.Nx), format='csr')
    M = Mlon + Mbas + Mlat

    rhs = -driving_k
    rhs[0] = 0.
    rhs[-1] = calving_ux # 0.5*self.delta*Hk[-1]**2
    
    M = M.asformat('lil') # seems to be happier about inserting entries in lil format
    M[0,0] = 1.
    M[0,1] = 0
    

    M[-1,-2] = self.Dxu[-1,-2]
    M[-1,-1] = self.Dxu[-1,-1]
    M = M.asformat('csr')
    
    norm = numpy.amax(numpy.abs(driving_k))
#    print "norm",norm
    self.resTau_k /= norm
    #print 'tau res k:', numpy.amax(numpy.abs(resTau_k))
       
    #resTau_k2 = (M.dot(uk) - rhs)/norm
   
    ukp1 = scipy.sparse.linalg.spsolve(M,rhs)
    
    self.resTau_solver = (M.dot(ukp1) - rhs)/norm
    self.noiseFloor = numpy.amax(numpy.abs(self.resTau_solver))
    #print 'tau solver res (noise floor):', self.noiseFloor

    self.u = ukp1
    self.resStress = numpy.amax(numpy.abs(self.resTau_k))
    
    
    
  def solveCont(self):
    Hk = self.H
    ukp1 = self.u
    
    a = self.a # non-dimensional accumulation
    #Hxk = self.DxH.dot(Hk)
    
   
    # build the up-wind matrix

    diagonals = []
    diag0 = ukp1*(ukp1 >= 0.)
    diagonals.append(diag0)
    diag1 = ukp1[0:-1]*(ukp1[0:-1] < 0.)
    diagonals.append(diag1)
    #print diagonals
#    print "diagonals,",diagonals    
    
    Mup = scipy.sparse.diags(diagonals, [0,1], shape=(self.Nx,self.Nx), format='csr')

    Mdt = scipy.sparse.diags([1./self.dt*numpy.ones((self.Nx),float)], [0],
                             shape=(self.Nx,self.Nx), format='csr')
    
    #print 'Mup:', Mup
    
    Hu = Hk*ukp1
    Hu2 = Hk[1:]*ukp1[0:-1]
    mask = ukp1[0:-1] < 0.
    Hu[0:-1][mask] = Hu2[mask]  #never happens in practice: u >= 0 everywhere
    
    M = self.Dxu.dot(Mup) + Mdt
    
    M = M.asformat('lil') # seems to be happier about inserting entries in lil format
    M[0,0] = self.DxH[0,0]
    M[0,1] = self.DxH[0,1]
    M = M.asformat('csr')
    

    rhs = a + self.HPrev/self.dt 
    rhs[0] = 0.
      
    Hkp1 = scipy.sparse.linalg.spsolve(M,rhs)
         
    #Hu = Hkp1*ukp1
    #resCont_kp1 = (Hkp1-self.HPrev)/self.dt + self.Dxu.dot(Hu) - a
    #resCont_kp1[0] = Hxk[0]
    
    self.H = Hkp1
    
    
    
  def step(self,outer):  
    for iterIndex in range(self.maxPicardIter):

      self.iterateOnViscosity()
      self.solveCont()
      self.updateXgAndFlotationMasks()
      
      self.plotIter(outer,iterIndex)
        
      #print resTau, tolerance
      #print "resStress", self.resStress
      #test residuals for convergence...
      if(self.resStress < self.toleranceInner):
        print 'Picard converged in %i iterations.'%(iterIndex+1)
        print 'resStress =', self.resStress, 'tol =',self.toleranceInner
        self.innerConverged = True
        return

      
    self.innerConverged = False
    print 'Warning:Picard iteration did not converge in', self.maxPicardIter, 'iterations. resStress:', self.resStress, 'tol:', self.toleranceInner

    
  def runTimeStepping(self,maxSteps,stepsPerWrite):

    converged = False
    
    self.updateXgAndFlotationMasks()
      
    for outer in range(maxSteps):
    
      self.uPrev = self.u
      self.HPrev = self.H
      self.xgPrev = self.xg
  
      self.step(outer)
  
      dH_dt = numpy.amax(numpy.abs(self.H-self.HPrev))/self.dt
      dxg_dt = numpy.abs(self.xg-self.xgPrev)/self.dt
  
      cfl = numpy.amax(numpy.abs(self.u/self.deltaX))*self.dt
      print "dH_dt = ",dH_dt, "dxg_dt = ",dxg_dt, "CFL = ",cfl
      print "xg = ", self.xg, "lambda_g = ", self.lambda_g
  
      if numpy.isnan(cfl):
        print "blew up!"
        exit(1)
      
#    print "time: ", solver.time+solver.dt, "iter: ", inner, "cfl:", cfl, "diffs: ", diffH, diffXg, diffU, solver.resStress
        
      if not self.innerConverged:
        print "Error: inner loop did not converge after %i steps!"%self.maxPicardIter
        print "Try reducing the goal CFL number."
#    print "diffU=,",diffU
        exit(1)
  

      self.time += self.dt
      print "outer:", outer, "time:", self.time
      if(dH_dt < self.tol_out_dH_dt and dxg_dt < self.tol_out_dxg_dt):
        converged = True
        break

  
      if(numpy.mod(outer,stepsPerWrite) == stepsPerWrite-1):
        self.writeResults()
  
    if(converged):
      print "The run converged."
      print "xg =",self.xg
    else:
      print "The run did not converge after %i steps."%maxSteps
      print "xg =",self.xg
  
    self.writeResults()

    if(not converged):
      exit(2)

