#!/usr/bin/python
import numpy
import Chebyshev
import numpy.linalg as linalg
  
class SteadyStateSolver:
  def __init__(self, Nx):
    self.Nx = Nx
    
    # Chebyshev operators
    self.cheb = Chebyshev.Chebyshev(Nx)
    self.intX_H = self.cheb.intX.copy()
    for index in range(Nx):
      self.intX_H[index,:] -= self.intX_H[-1,:]
      
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
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    self.n = 3. # Glen's flow parameter

    self.lambda_0 = 2   # wavelength of bedrock bump (m)
    self.m_0  = .5  # maximum bed obstacle slope (no unit)
    
    self.eps_s = 1e-3
    
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
    
    self.delta = 1.0 - self.rho_i/self.rho_w
        
  # call each time H and xg are updated to compute b, Hf, Hx, bx and Np
  def updateH(self, H, xg):
    Dx = self.cheb.Dx/xg
    x = xg*self.cheb.x

    self.H = H
    self.xg = xg

    (self.b,self.bx) = self.computeB(x,self)
    self.Hf = self.maxS(self.b)/(1-self.delta)

    self.Hx = numpy.dot(Dx,self.H)
    #self.bx = numpy.dot(Dx,self.b)
    self.Np = self.H*(self.maxS(1.0 - self.Hf/self.H))**self.p
  
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
    f = self.Hf - newH
    fx_xg = numpy.dot(Dx[-1,:],f)
    newXg = xg - f[-1]/fx_xg
    
    return newXg


  def iterateImplicitTimeStep(self):
    # u Hx + ux H = a
    xg = self.xg
    sigma = self.cheb.x
    x = xg*self.cheb.x
    Dx = self.cheb.Dx/xg
    
    dxg_dt = (xg - self.oldXg)/self.dt
    
    M = numpy.zeros((self.Nx+1,self.Nx+1))
    rhs = numpy.zeros((self.Nx+1,1))
    
    # dH/dt - sigma dxg_dt Hx + u Hx + ux H = a
    M =  numpy.diag(1./self.dt  + self.ux) + numpy.dot(numpy.diag(self.u - sigma*dxg_dt),Dx)
    rhs = self.a + self.oldH/self.dt
    # no boundary conditions should be necessary -- ux H = a at the ice divide
    newH = linalg.solve(M,rhs)
    
    newHx = numpy.dot(Dx,newH)
    if(self.plot):
      import matplotlib.pyplot as plt
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,5)
      ax.cla()
      plt.plot(x,newH-self.oldH, 'b')
      ax = plt.subplot(2,3,6)
      ax.cla()
      plt.plot(x,(newH-self.oldH)/self.dt, 'b', x,-sigma*dxg_dt*newHx, 'r', 
               x, self.u*newHx, 'k', x, self.ux*newH, 'g')
#      if(self.plotContinuous):
#        plt.draw()
#      else:
#        plt.show()
    return newH
    
#  def iterateTowardSteadyState(self,alpha):
#    # u Hx + ux H = a
#    xg = self.xg
#    x = xg*self.cheb.x
#    Dx = self.cheb.Dx/xg
#        
#    M = numpy.zeros((self.Nx+1,self.Nx+1))
#    rhs = numpy.zeros((self.Nx+1,1))
#    
#    # dH/dt - sigma dxg_dt Hx + u Hx + ux H = a
#    M =  numpy.diag(self.ux) + numpy.dot(numpy.diag(self.u),Dx)
#    rhs = self.a*numpy.ones(self.H.shape)
#    # no boundary conditions should be necessary -- ux H = a at the ice divide
#    steadyH = linalg.solve(M,rhs)
#    
#    #scale = numpy.amax(self.oldH)/numpy.amax(steadyH)
#    newH = (1-alpha)*self.oldH + alpha*steadyH
#    
#    newHx = numpy.dot(Dx,newH)
#    if(self.plot):
#      plt.subplot(2,3,5)
#      plt.plot(x,newH-self.oldH, 'b')
#      plt.title('xg=%f'%xg)
#      plt.subplot(2,3,6)
#      plt.plot(x, self.u*newHx, 'r', 
#               x, self.ux*newH, 'b', 
#               x, self.a*numpy.ones(self.H.shape), 'g', 
#               x, self.u*newHx + self.ux*newH - self.a, 'k')
#    return newH

  def iterateOnViscosity(self):
    # make local copies of common variables for convenience
    xg = self.xg
    x = xg*self.cheb.x
    Dx = self.cheb.Dx/xg
    n = self.n
    
    # these variables should be computed by updateH()
    # before iteration begins
    H = self.H
    Hx = self.Hx
    b = self.b
    bx = self.bx
    Hf = self.Hf
    Np = self.Np
    
    uk = self.u # previous iteration
    uxk = self.ux


    #tau_l = (H nu ux)_x
    #nu = 4*epsilon*|ux|^(1/n - 1)
    abs_ux = self.absS(uxk)
    if(self.useLongi):
      # tauLCoeff = H nu
      tauLCoeff = 4.*self.epsilon*(abs_ux)**(1./n-1.)*H
    else:
      tauLCoeff = numpy.zeros(H.shape)

    #tau_b = -gamma*(N^n/(Kappa*|u| + N^n))^(1/n) u/|u|
    abs_u = self.absS(uk)
    if(self.useSchoofBasal):
      tauBCoeff = -self.gamma/abs_u
    else:
      tauBCoeff = -self.gamma*Np*(abs_u/(self.Kappa*abs_u + Np**n))**(1./n)/abs_u


    #tau_d = -H (H-b)_x
    tauD = -H*(Hx-bx)

    tauLk = numpy.dot(Dx,tauLCoeff*uxk)
    tauBk = tauBCoeff*uk
    
    residual = tauLk + tauBk + tauD
    
    M = numpy.dot(Dx,numpy.dot(numpy.diag(tauLCoeff),Dx)) + numpy.diag(tauBCoeff)
    rhs = -tauD
    
    # velocity is zero at the ice divide
    M[0,:] = 0.
    M[0,0] = 1.
    rhs[0] =  0.
    
    # H nu ux = 0.5*delta*H^2 at x = xg
    tauLCoeff_xg = 4.*self.epsilon*(abs_ux[-1])**(1./n-1.)*H[-1]
    M[-1,:] = Dx[-1,:]
    rhs[-1] = 0.5*self.delta*H[-1]**2/tauLCoeff_xg
    
    
    ukp1 = linalg.solve(M,rhs)
    uxkp1 = numpy.dot(Dx,ukp1)
    
    solver_res = numpy.dot(M,ukp1) - rhs
    residual_check = numpy.dot(M,uk) - rhs
    self.noiseFloor = numpy.max(numpy.abs(solver_res[1:-1]))

    if(self.plot):
      import matplotlib.pyplot as plt
      plt.figure(1, figsize=(16, 12))
      ax = plt.subplot(2,3,1)
      ax.cla()
      oldX = self.cheb.x*self.oldXg
      (oldB,oldBx) = self.computeB(oldX,self)
      plt.plot(x,H-b, 'r', x, -b, 'k', x, Hf-b, 'g',oldX,self.oldH-oldB,'b')
      ax = plt.subplot(2,3,2)
      ax.cla()
      plt.plot(x[1:-1],solver_res[1:-1],'b', x[1:-1], residual[1:-1], 'r',
               x[1:-1], residual_check[1:-1], 'k--')
      ax = plt.subplot(2,3,3)
      ax.cla()
      plt.plot(x,tauLk, 'b', x, tauBk, 'r', x, tauD, 'g', x, residual, 'k')
      ax = plt.subplot(2,3,4)
      ax.cla()
      plt.plot(x,uk, 'b', x, ukp1, 'r')
      #plt.subplot(2,3,5)
      #plt.plot(x,uxk,'b',x,uxkp1,'r',x,rhs[-1]*numpy.ones(x.shape),'g')
#      if(self.plotContinuous):
#        plt.draw()
#      else:
#        plt.show()
      
      
    self.resStress = numpy.max(numpy.abs(residual[1:-1]))
    
    self.u = ukp1
    self.ux = uxkp1
    
