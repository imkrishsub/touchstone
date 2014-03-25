# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 15:37:32 2013

@author: xylar
"""
from pprint import pprint
import numpy
import scipy.sparse
import scipy.sparse.linalg
import numpy.linalg

import matplotlib.pyplot as plt

class Solver:
  def __init__(self, Nx,deltaX,xc,p,A,linearBed,linearSlope,useGLP):
    if(linearBed):
      self.computeB = self.computeBLinear
    else:
      self.computeB = self.computeBPoly
    self.Nx = Nx
    self.deltaX = deltaX
    self.xc = xc
    self.xu = numpy.linspace(0,xc,Nx)
    self.xH = self.xu-0.5*self.deltaX
    self.rho_i = 900.0 # kg/m^3 ice density
    self.rho_w = 1000.0 # kg/m^3 water density
    self.delta = 1.0 - self.rho_i/self.rho_w
    self.C = numpy.array(7.624e6)
    self.W = 10.
    self.a = numpy.array(1.0)
    self.p = numpy.array(p)
    self.A = numpy.array(A)
    self.g = 9.8 # m/s^2 gravity acceleration
    self.n = 3. # Glen's flow parameter
    self.linearSlope = linearSlope

    self.Ab = 3.1688e-24 #Pa^-3 s^-1
    self.lambda_0 = 2   # wavelength of bedrock bump (m)
    self.m_0  = .5  # maximum bed obstacle slope (no unit)

    self.eps_s = 1e-10
    self.tolerance = 1e-10

    self.plotContinuous = True
    self.plot = True
    self.useLongi = True
    self.useChannel = False
    self.useSchoofBasal = False
    self.useGLP = False
    self.transient = True
    self.maxPicardIter = 200
    self.maxPicardInner = 1
    
    sPerY = 365.25*24.*3600. # number of seconds per year
    self.aBar = .3/sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.WBar = 10000. # m channel width scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor    
    self.tBar = self.xBar/self.uBar    
    
    self.NBar = self.rho_i*self.g*self.HBar
    self.taudBar = self.NBar*self.HBar/self.xBar
    self.taubBar = self.C*self.uBar**(1./self.n)
    self.Kappa = (self.m_0*self.uBar)/(self.lambda_0*self.Ab)/self.NBar**self.n
    self.gamma = self.taubBar/self.taudBar

    self.taulBar = (A**(-1./self.n))*(self.uBar/self.xBar)**(1./self.n)*self.HBar/self.xBar
    self.epsilon = self.taulBar/(2*self.taudBar)

    self.tauWbar = self.HBar/self.WBar*(5*self.uBar/(2*A*self.WBar))**(1./self.n)
    self.omega = self.tauWbar/self.taudBar
    
    self.timeCentering = 1.
    
    print 'gamma, Kappa, epsilon, omega:', self.gamma, self.Kappa, self.epsilon, self.omega
    
#    print "self:"
#    pprint(self.__dict__)    
    
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
    print 'b:', self.b
    self.Hf = self.maxS(self.b)/(1-self.delta)
    self.bx = self.DxH.dot(self.b)


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

  def computeSx(self,H,floatingMaskH,floatingMaskU,glIndices,lambda_g):
    s = H-self.b
    s_floating = self.delta*H
    s[floatingMaskH] = s_floating[floatingMaskH]
    
    sx = self.DxH.dot(s)
    sx_floating = self.DxH.dot(s_floating)
    sx[floatingMaskU] = sx_floating[floatingMaskU]
    
    if(self.useGLP):
        sx[glIndices] = self.lambda_g*sx[glIndices] \
        + (1. - self.lambda_g)*sx_floating[glIndices]
    
    return (s,sx)
    
    
    
  def getFlotationMasks(self, H):
    fPattyn = self.Hf/H
    
    groundedMaskH = fPattyn < 1.
    floatingMaskH = fPattyn >= 1.
    
    groundedMaskU = numpy.zeros(groundedMaskH.shape,bool)
    groundedMaskU[0:-1] = numpy.logical_and(groundedMaskH[1:],groundedMaskH[0:-1])
    
    floatingMaskU = numpy.ones(floatingMaskH.shape,bool)
    floatingMaskU[0:-1] = numpy.logical_and(floatingMaskH[1:],floatingMaskH[0:-1])
    
    glIndices = numpy.nonzero(numpy.logical_xor(groundedMaskH[1:],groundedMaskH[0:-1]))[0]
    #print 'gl at index (or indices):', glIndices
        
    lambda_g = (1. - fPattyn[glIndices])/(fPattyn[glIndices+1] - fPattyn[glIndices])
    
    
    return (groundedMaskH,floatingMaskH,groundedMaskU,floatingMaskU,glIndices,lambda_g)
 

  def iteratePicardTau(self, Hk, uk):
    p = self.p
    n = self.n # Glen's flow parameter
    #xH = self.xH
    xu = self.xu
    Hf = self.Hf    
    
    uxk = self.Dxu.dot(uk)
    abs_ux = self.absS(uxk)
    nuk = self.epsilon*(abs_ux)**(1./n-1.)
    
    fPattyn = Hf/Hk
    
    groundedMaskH = fPattyn < 1.
    floatingMaskH = fPattyn >= 1.
  
    glIndices = numpy.nonzero(numpy.logical_xor(groundedMaskH[1:],groundedMaskH[0:-1]))[0]
        
    lambda_g = (1. - fPattyn[glIndices])/(fPattyn[glIndices+1] - fPattyn[glIndices])
  
    groundedMaskU = numpy.zeros(groundedMaskH.shape,bool)
    groundedMaskU[0:-1] = numpy.logical_and(groundedMaskH[1:],groundedMaskH[0:-1])
    
    floatingMaskU = numpy.ones(floatingMaskH.shape,bool)
    floatingMaskU[0:-1] = numpy.logical_and(floatingMaskH[1:],floatingMaskH[0:-1])
    
  
    # setting up tau_b  
     
    abs_u = self.absS(uk)
    if(self.useSchoofBasal):
      basalCoeff = -self.gamma*abs_u**(1./n-1)
    else:
      Nk = Hk*(self.maxS(1.0 - self.Hf/Hk))**p
      self.N = Nk
      Nk_u = self.AvgH.dot(Nk)
      Nk_u[floatingMaskU] = 0.
      #Nk_u[glIndices] = numpy.maximum(Nk[glIndices],Nk[glIndices+1])
      self.Nu = Nk_u
      gammaEff = self.gamma*(Nk_u**n/(self.Kappa*abs_u + Nk_u**n))**(1./n)
      self.basalu = -self.gamma*abs_u**(1./n-1)*uk
      self.basalN = -self.gamma/(self.Kappa)**(1/n)*Nk_u*abs_u**(-1)*uk
      basalCoeff = -gammaEff*abs_u**(1./n-1)
                 
    #print self.gamma
    if(self.useGLP):
        basalCoeff[glIndices] *= lambda_g #*2/(1+p)
   
    basalCoeff[floatingMaskU] = 0.
  
  
    # setting up Tau_l        
    
    longiCoeff = 4*Hk*nuk
    if(not self.useLongi):
      longiCoeff[self.groundedMaskU] = 0.
      

    # setting up Tau_w

    if(not self.useChannel):
        lateralEff = 0.
    else:
        lateralEff = self.omega*Hk*self.W**(1./n+1)
      
    lateralCoeff = -lateralEff*abs_u**(1./n-1)

  
    # setting up Tau_d
  
    (sk,sxk) = self.computeSx(Hk,floatingMaskH,floatingMaskU,glIndices,lambda_g)      
    Hk_u = self.AvgH.dot(Hk)
    driving_k = -Hk_u*sxk
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

    resTau_k = longi_k + basal_k + driving_k + lateral_k
    resTau_k[0] = uk[0]
    resTau_k[-1] = uxk[-1] - calving_ux
    
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
    resTau_k /= norm
    #print 'tau res k:', numpy.amax(numpy.abs(resTau_k))
       
    resTau_k2 = (M.dot(uk) - rhs)/norm
   
    ukp1 = scipy.sparse.linalg.spsolve(M,rhs)
    
    resTau_kp1 = (M.dot(ukp1) - rhs)/norm
    self.noiseFloor = numpy.amax(numpy.abs(resTau_kp1))
    #print 'tau solver res (noise floor):', self.noiseFloor
    
    if(self.plot):
      fig = plt.figure(1)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xu,uk, 'b', xu,ukp1, 'r')        
      fig = plt.figure(2)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xu,longi_k, 'b', xu, basal_k, 'r', xu, driving_k, 'g', xu, resTau_k*norm, 'k',
               xu, resTau_k2*norm, 'k--')

      if(self.plotContinuous):
        plt.draw()

    self.u = ukp1
    self.resStress = numpy.linalg.norm(resTau_kp1)
#    self.resStress = numpy.max(numpy.abs(resTau_k2))
    #return (ukp1,numpy.max(numpy.abs(resTau_k2)))
    
    
    
  def solveCont(self, Hk, ukp1, HPrev, deltaT, firstIter):
    a = self.a # non-dimensional accumulation
    Hxk = self.DxH.dot(Hk)
    
    if(not self.transient):
      self.timeCentering = 1.0
   
    # build the up-wind matrix

    diagonals = []
    diag0 = ukp1*(ukp1 >= 0.)
    diagonals.append(diag0)
    diag1 = ukp1[0:-1]*(ukp1[0:-1] < 0.)
    diagonals.append(diag1)
    #print diagonals
#    print "diagonals,",diagonals    
    
    Mup = scipy.sparse.diags(diagonals, [0,1], shape=(self.Nx,self.Nx), format='csr')

    Mdt = scipy.sparse.diags([1./deltaT*numpy.ones((self.Nx),float)], [0],
                             shape=(self.Nx,self.Nx), format='csr')
    
    #print 'Mup:', Mup
    
    Hu = Hk*ukp1
    Hu2 = Hk[1:]*ukp1[0:-1]
    mask = ukp1[0:-1] < 0.
    Hu[0:-1][mask] = Hu2[mask]  #never happens in practice: u >= 0 everywhere
    
    if(firstIter):
      self.prevHux = self.Dxu.dot(Hu)
    
    M = self.timeCentering*self.Dxu.dot(Mup) + float(self.transient)*Mdt
    
    M = M.asformat('lil') # seems to be happier about inserting entries in lil format
    M[0,0] = self.DxH[0,0]
    M[0,1] = self.DxH[0,1]
    M = M.asformat('csr')
    

#    rhs = a + self.transient*HPrev/deltaT - (1-self.timeCentering)*self.prevHux
    rhs = a + self.transient*HPrev/deltaT 
    rhs[0] = 0.
      
    Hkp1 = scipy.sparse.linalg.spsolve(M,rhs)
         
    resCont_kp1 = float(self.transient)*(Hkp1-HPrev)/deltaT + self.Dxu.dot(Hu) - a
    resCont_kp1[0] = Hxk[0]
    
    return Hkp1
    
    
    
  def step(self, HGuess, uGuess, HPrev, deltaT):
    if(self.plotContinuous):
      plt.ion()
  
    Hk = HGuess
    uk = uGuess
#    (self.groundedMaskH,self.floatingMaskH,self.groundedMaskU,self.floatingMaskU,
#     self.glIndices,self.lambda_g) \
#      = self.getFlotationMasks(Hk)
#
#    (groundedMaskHPrev,floatingMaskHPrev,groundedMaskUPrev,
#     floatingMaskUPrev,
#     self.glIndicesPrev,self.lambda_gPrev) \
#      = self.getFlotationMasks(HPrev)
#    (self.sPrev,sxPrev) = self.computeSx(HPrev,floatingMaskHPrev,
#        floatingMaskUPrev,self.glIndicesPrev,self.lambda_gPrev)
   
    for iterIndex in range(self.maxPicardIter):
#      (groundedMaskH,floatingMaskH,groundedMaskU,floatingMaskU,glIndices,lambda_g) \
#        = self.getFlotationMasks(Hk)

      for inner in range(self.maxPicardInner):
        self.iteratePicardTau(Hk, uk)
        uk = self.u # it's really ukp1
        
#      if (iterIndex==1):
#          firstIter = True
#      else:
#          firstIter = False
      HPrev = Hk
      Hkp1 = self.solveCont(Hk, uk, HPrev, deltaT, (iterIndex == 0))
        
#      tolerance = numpy.maximum(self.tolerance,2*self.noiseFloor)
      tolerance = 1e-8  
      #print 'resTau, tol:', resTau, tolerance
      if(self.plot and not self.plotContinuous):
        plt.ioff()
        plt.show()  
        
      #print resTau, tolerance
      #test residuals for convergence...
      if((self.resStress < tolerance) and iterIndex > 0):
        print 'Picard converged in %i iterations.'%iterIndex
        print 'resStress, tol:', self.resStress, tolerance
        self.innerConverged = True
        if(self.plot):
          plt.ioff()
          plt.show()
        return (Hkp1,uk)

      Hk = Hkp1
      #uk = ukp1
      
    print 'Warning:Picard iteration did not converge in', self.maxPicardIter, 'iterations. resStress:', self.resStress, 'tol:', tolerance
    #raise ValueError('Picard did not converge in %i iterations.'%self.maxPicardIter)
    return (Hkp1,uk)
    
    
    
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
        if(deltaH < self.tolerance):
          break
      #print "deltaH:", Hi-HPrev, Hi, Hip1
      Hip1 = Hi
      H[xIndex] = Hi
      u[xIndex] = umid
    
  
    return (H,u)


    
  def updateXg(self,H,Hf):
    # this function updates xg at each new time step, i.e when newH is updated
     diff = Hf-H
 
     (self.groundedMaskH,self.floatingMaskH,self.groundedMaskU,self.floatingMaskU,
      self.glIndices,self.lambda_g) \
       = self.getFlotationMasks(H)    
        
     s = diff[self.glIndices]/(diff[self.glIndices]-diff[self.glIndices-1])
     xH = self.xH
     Xg = s*xH[self.glIndices-1]+(1-s)*xH[self.glIndices]
        
     return Xg 
        
      


