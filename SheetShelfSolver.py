#!/usr/bin/python
import numpy
import Chebyshev
import numpy.linalg as linalg
  
class SheetShelfSolver:
  def __init__(self, Nx):
    self.Nx = Nx
    
    # Chebyshev operators
    self.cheb = Chebyshev.Chebyshev(Nx)
      
    # default values for constants
    self.rho_i = 900.0 # kg/m^3 ice density
    self.rho_w = 1000.0 # kg/m^3 water density
    self.C = 7.624e6
    self.a = 1.0
    self.plotContinuous = True
    self.plot = True
    self.useLongi = True
    self.useSchoofBasal = False
    self.Ab = 3.1688e-24 #Pa^-3 s^-1
    self.A = 1e-25 #Pa^-3 s^-1, a typical value
    
    self.p = 0.0
    
    sPerY = 365.25*24.*3600. # number of seconds per year
    self.g = 9.8 # m/s^2 gravity acceleration
    self.aBar = .3/sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m height scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    self.n = 3. # Glen's flow parameter
    self.WBar = 1000. # m width scaling factor

    self.lambda_0 = 2   # wavelength of bedrock bump (m)
    self.m_0  = .5  # maximum bed obstacle slope (no unit)
    
    self.eps_s = 1e-3
    
    self.xc = 2.4
    self.useChannelWidth = False
    
    self.computeNondimConstants()
    
  # call this function every time a constant (e.g., A, C, n, Ab, m, lambda_0, m_0, etc.)
  # is updated or changed from the default values,
  def computeNondimConstants(self):
    n = self.n
    NBar = self.rho_i*self.g*self.HBar
    tauLBar = self.A**(-1./n)*self.HBar*self.uBar**(1./n)/self.xBar**(1./n + 1)
    tauBBar = self.C*self.uBar**(1./n)
    tauDBar = self.rho_i*self.g*self.HBar**2/self.xBar
    self.Kappa = self.m_0/(self.lambda_0*self.Ab)*self.uBar/NBar**n
    self.epsilon = tauLBar/(2*tauDBar)
    self.gamma = tauBBar/tauDBar
    
    tauWBar = self.HBar/self.WBar*(5.*self.uBar/(2.*self.A*self.WBar))**(1./n)
    self.omega = tauWBar/tauDBar
    
    self.delta = 1.0 - self.rho_i/self.rho_w
        
  # call each time H and xg are updated to compute b, Hf, Hx, bx and Np
  def updateH(self, H, xg):
    Nx = self.Nx
    xc = self.xc
    if(xg >= xc):
      print "Error: xg=%f >= xc=%f. Exiting..."%(xg,xc)
      exit(1)
    DxSheet = self.cheb.Dx/xg
    xSheet = xg*self.cheb.x
    DxShelf = self.cheb.Dx/(xc-xg)
    xShelf = xg + (xc-xg)*self.cheb.x
    
    x = numpy.zeros((2*Nx))
    x[0:Nx] = xSheet
    x[Nx:2*Nx] = xShelf

    self.H = H
    self.xg = xg
    self.xc = xc
    self.x = x

    HxSheet = numpy.dot(DxSheet,self.H[0:Nx])
    HxShelf = numpy.dot(DxShelf,self.H[Nx:2*Nx])

    (self.b,bx) = self.computeB(x,self)
    self.Hf = self.maxS(self.b[0:Nx])/(1-self.delta)
    self.Np = self.H[0:Nx]*(self.maxS(1.0 - self.Hf/self.H[0:Nx]))**self.p

    self.s = numpy.zeros(x.shape)
    self.s[0:Nx] = self.H[0:Nx]-self.b[0:Nx]
    self.s[Nx:2*Nx] = self.delta*self.H[Nx:2*Nx]

    self.sx = numpy.zeros(x.shape)
    self.sx[0:Nx] = HxSheet-bx[0:Nx]
    self.sx[Nx:2*Nx] = self.delta*HxShelf
    
    if(self.useChannelWidth):
      (self.W,self.Wx) = self.computeW(x,self)
  
  def absS(self,x):
    return numpy.sqrt(x**2 + self.eps_s**2)
    
  def absSx(self,x):
    return x/self.absS(x)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)
    
  def maxSx(self,x):
    return 0.5*(self.absSx(x) + 1.0)
    
  def newtonStepXg(self, newH):
    xg = self.xg
    Dx = self.cheb.Dx/xg
    
    # Newton's method to find the new xg
    # xg^(k+1) = xg^k - f(xg^k)/f'(xg^k)
    # f = Hf-H
    f = self.Hf - newH[0:self.Nx]
    fx_xg = numpy.dot(Dx[-1,:],f)
    newXg = xg - f[-1]/fx_xg
    
    return newXg


  def iterateImplicitTimeStep(self):
    # u Hx + ux H = a
    xg = self.xg
    sigma = self.cheb.x
    Nx = self.Nx
    xc = self.xc
    DxSheet = self.cheb.Dx/xg
    DxShelf = self.cheb.Dx/(xc-xg)
    
    x = self.x
    
    Dx = numpy.zeros((2*Nx,2*Nx))
    Dx[0:Nx,0:Nx] = DxSheet
    Dx[Nx:2*Nx,Nx:2*Nx] = DxShelf

    dxg_dt = (xg - self.oldXg)/self.dt    
    
    movingGridTerm = numpy.zeros(x.shape)
    movingGridTerm[0:Nx] = -sigma*dxg_dt
    movingGridTerm[Nx:2*Nx] = -(1.-sigma)*dxg_dt
    
    diag1 = 1./self.dt  + self.ux
    if(self.useChannelWidth):
      diag1 += (self.u*self.Wx)/self.W

    # dH/dt - sigma dxg_dt Hx + u Hx + ux H + u H Wx/W = a (in sheet)
    M =  numpy.diag(diag1) + numpy.dot(numpy.diag(self.u + movingGridTerm),Dx)
    rhs = self.a + self.oldH/self.dt
    
    # have to replace a row to get continuity in H at xg
    M[Nx,:] = 0.
    M[Nx,Nx-1] = 1.
    M[Nx,Nx] = -1.
    rhs[Nx] = 0.
    
    # no boundary conditions should be necessary -- ux H = a at the ice divide
    newH = linalg.solve(M,rhs)
    
    if(self.plot):
      newHx = numpy.dot(Dx,newH)
      import matplotlib.pyplot as plt
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,6)
      ax.cla()
      plt.plot(x,(newH-self.oldH)/self.dt, 'b', x,movingGridTerm*newHx, 'r', 
               x, self.u*newHx, 'k', x, self.ux*newH, 'g')
    return newH

  def iterateOnViscosity(self):
    # make local copies of common variables for convenience
    xg = self.xg
    Nx = self.Nx
    xc = self.xc
    DxSheet = self.cheb.Dx/xg
    DxShelf = self.cheb.Dx/(xc-xg)
    
    x = self.x
    
    Dx = numpy.zeros((2*Nx,2*Nx))
    Dx[0:Nx,0:Nx] = DxSheet
    Dx[Nx:2*Nx,Nx:2*Nx] = DxShelf

    n = self.n
    
    # these variables should be computed by updateH()
    # before iteration begins
    H = self.H
    #b = self.b
    Hf = self.Hf
    Np = self.Np
    
    uk = self.u # previous iteration
    uxk = self.ux


    #tau_l = (H nu ux)_x
    #nu = 4*epsilon*|ux|^(1/n - 1)
    abs_ux = self.absS(uxk)
    tauLCoeff = 4.*self.epsilon*(abs_ux)**(1./n-1.)*H
    if(not self.useLongi):
      tauLCoeff[0:Nx] = 0.

    #tau_b = -gamma*(N^n/(Kappa*|u| + N^n))^(1/n) u/|u|
    abs_u = self.absS(uk)
    tauBCoeff = numpy.zeros(x.shape)
    if(self.useSchoofBasal):
      tauBCoeff[0:Nx] = -self.gamma*abs_u[0:Nx]**(1/n-1)
    else:
      tauBCoeff[0:Nx] = -self.gamma*Np*(abs_u[0:Nx]/(self.Kappa*abs_u[0:Nx] + Np**n))**(1./n)/abs_u[0:Nx]

    #tau_w = -omega*H*W^(-1-1/n)*|u|^(1/n-1) u
    if(self.useChannelWidth):
      tauWCoeff = -self.omega*H*self.W**(-1.-1./n)*abs_u**(1/n-1)
    else:
      tauWCoeff = numpy.zeros(x.shape)

    #tau_d = -H sx
    tauD = -H*self.sx

    tauLk = numpy.dot(Dx,tauLCoeff*uxk)
    tauBk = tauBCoeff*uk
    tauWk = tauWCoeff*uk
    
    residual = tauLk + tauBk + tauD + tauWk
    
    M = numpy.dot(Dx,numpy.dot(numpy.diag(tauLCoeff),Dx)) + numpy.diag(tauBCoeff + tauWCoeff)
    rhs = -tauD
    
    # velocity is zero at the ice divide
    M[0,:] = 0.
    M[0,0] = 1.
    rhs[0] =  0.
    
    # H nu ux = 0.5*delta*H^2 at x = xc
    tauLCoeff_xg = 4.*self.epsilon*(abs_ux[-1])**(1./n-1.)*H[-1]
    M[-1,:] = Dx[-1,:]
    rhs[-1] = 0.5*self.delta*H[-1]**2/tauLCoeff_xg
    
    # have to replace rows to get continuity in u and ux at xg
    M[Nx-1,:] = 0.
    M[Nx-1,Nx-1] = 1.
    M[Nx-1,Nx] = -1.
    rhs[Nx-1] = 0.
    M[Nx,0:Nx] = DxSheet[Nx-1,:]
    M[Nx,Nx:2*Nx] = -DxShelf[0,:]
    rhs[Nx] = 0.
        
    ukp1 = linalg.solve(M,rhs)
    uxkp1 = numpy.dot(Dx,ukp1)
    
    #self.noiseFloor = numpy.max(numpy.abs(solver_res[1:-1]))

    if(self.plot):
      import matplotlib.pyplot as plt
      residual_check = numpy.dot(M,uk) - rhs
      solver_res = numpy.dot(M,ukp1) - rhs
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,1)
      ax.cla()
      oldX = self.cheb.x*self.oldXg
      (oldB,oldBx) = self.computeB(oldX,self)
      plt.plot(x,self.s, 'r', x, -self.b, 'k', x[0:Nx], Hf-self.b[0:Nx], 'g',oldX,self.oldH[0:Nx]-oldB,'b')
      ax = plt.subplot(2,3,2)
      ax.cla()
      plt.plot(x[1:-1],solver_res[1:-1],'b', x[1:-1], residual[1:-1], 'r',
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
    
