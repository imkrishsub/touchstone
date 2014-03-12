# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 15:37:32 2013

@author: xylar
"""

import numpy
import scipy.sparse
import scipy.sparse.linalg
import numpy.linalg

import matplotlib.pyplot as plt

class Solver:
  def __init__(self, Nx,deltaX,p,A,linearBed):
    if(linearBed):
      self.computeB = self.computeBLinear
    else:
      self.computeB = self.computeBPoly
    self.Nx = Nx
    #self.Lx = Lx
    #self.deltaX = Lx/(Nx-1.5)
    self.deltaX = deltaX
    self.Lx = (Nx-1.5)*deltaX
    print "Lx =", self.Lx
    self.xu = self.deltaX*numpy.arange(Nx)
    self.xH = self.xu-0.5*self.deltaX
    self.rho_i = 900.0 # kg/m^3 ice density
    self.rho_w = 1000.0 # kg/m^3 water density
    self.delta = 1.0 - self.rho_i/self.rho_w
    self.C = numpy.array(7.624e6)
    self.a = numpy.array(1.0)
    self.p = numpy.array(p)
    self.A = numpy.array(A)
    self.g = 9.8 # m/s^2 gravity acceleration
    self.n = 3. # Glen's flow parameter

    self.Ab = 3.1688e-24 #Pa^-3 s^-1
    self.lambda_0 = 2   # wavelength of bedrock bump (m)
    self.m_0  = .5  # maximum bed obstacle slope (no unit)

    self.eps_s = 1e-10
    self.tolerance = 1e-10

    self.plotContinuous = True
    self.plot = True
    self.useLongi = True
    self.useSchoofBasal = False
    self.transient = True
    self.maxPicardIter = 200
    self.maxPicardInner = 1
    
    sPerY = 365.25*24.*3600. # number of seconds per year
    self.aBar = .3/sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    
    self.NBar = self.rho_i*self.g*self.HBar
    self.taudBar = self.NBar*self.HBar/self.xBar
    self.taubBar = self.C*self.uBar**(1./self.n)
    self.Kappa = (self.m_0*self.uBar)/(self.lambda_0*self.Ab)/self.NBar**self.n
    self.gamma = self.taubBar/self.taudBar

    self.taulBar = (A**(-1./self.n))*(self.uBar/self.xBar)**(1./self.n)*self.HBar/self.xBar
    self.epsilon = self.taulBar/(2*self.taudBar)
    
    self.timeCentering = 0.5
    
    print 'gamma, Kappa, epsilon:', self.gamma, self.Kappa, self.epsilon
    
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

    self.b = self.computeB(self.xH)
    self.Hf = self.maxS(self.b)/(1-self.delta)
    self.bx = self.DxH.dot(self.b)


   
  def computeBLinear(self, x):
    shoofx = 750000.
    bx = 778.5*self.xBar/shoofx/self.HBar
    b0 = -720./self.HBar
    eps_b = 1e-3
    absx = numpy.sqrt(x**2 + eps_b**2)
    b = b0 + bx*absx #bx*numpy.abs(x)
    return b
    
  def computeBPoly(self, x):
    schoofx = 750000.
    xs = x*self.xBar/schoofx
    b = -(729 - 2184.8*xs**2 + 1031.72*xs**4 - 151.72*xs**6)/self.HBar  
    return b


  def absS(self,x):
    return numpy.sqrt(x**2 + self.eps_s**2)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)

#  def dxH(self,H):
#    Hx = numpy.zeros(H.shape,float)
#    Hx[0:-1] = (H[1:]-H[0:-1])/self.deltaX
#    return Hx
#
#  def dxu(self,u):
#    ux = numpy.zeros(u.shape,float)
#    ux[1:] = (u[1:]-u[0:-1])/self.deltaX
#    return ux
#
#
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
    
    sx = self.DxH.dot(s)
    sx_floating = self.DxH.dot(s_floating)
    sx[floatingMaskU] = sx_floating[floatingMaskU]
    sx[glIndices] = self.lambda_g*sx[glIndices] \
      + (1. - self.lambda_g)*sx_floating[glIndices]
      
    s[floatingMaskH] = s_floating[floatingMaskH]
#    plt.figure(5)
#    plt.plot(self.xH,s, 'b', self.xH, s_floating, 'r', self.xH, H-self.b, 'm')
#    plt.figure(6)
#    plt.plot(self.xu,sx, 'b', self.xu, sx_floating, 'r', self.xu, self.DxH.dot(H-self.b), 'm')
#    plt.show()
    
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
    uxk = self.Dxu.dot(uk)
    abs_ux = self.absS(uxk)
    nuk = self.epsilon*(abs_ux)**(1./n-1.)
    
    floatingMaskU = self.floatingMaskU
    floatingMaskH = self.floatingMaskH
    glIndices = self.glIndices
    lambda_g = self.lambda_g
     
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
    
    basalCoeff[glIndices] *= lambda_g #*2/(1+p)
    basalCoeff[floatingMaskU] = 0.
  
    longiCoeff = 4*Hk*nuk
    if(not self.useLongi):
      longiCoeff[self.groundedMaskU] = 0.
  

    (sk,sxk) = self.computeSx(Hk,floatingMaskH,floatingMaskU,glIndices,lambda_g)
          
    Hk_u = self.AvgH.dot(Hk)
    
    driving_k = -Hk_u*sxk
#    if(hasattr(self,'driving')):
#      deltaDriving = driving_k-self.driving
#    else:
#      deltaDriving = numpy.zeros(driving_k.shape)
    self.driving = driving_k
    
    calving_ux = (self.delta/(8*self.epsilon)*Hk[-1])**n
    #print 'calving_ux:', calving_ux, uxk[-1]
    
    
    longi_k = self.DxH.dot(longiCoeff*uxk)
#    if(hasattr(self,'longi')):
#      deltaLongi = longi_k-self.longi
#    else:
#      deltaLongi = numpy.zeros(longi_k.shape)
    self.longi = longi_k
    basal_k = basalCoeff*uk
#    if(hasattr(self,'basal')):
#      deltaBasal = basal_k-self.basal
#    else:
#      deltaBasal = numpy.zeros(basal_k.shape)
    self.basal = basal_k
    resTau_k = longi_k + basal_k + driving_k
    resTau_k[0] = uk[0]
    resTau_k[-1] = uxk[-1] - calving_ux

    Mlon = scipy.sparse.diags([longiCoeff],[0],shape=(self.Nx,self.Nx), format='csr')
    Mlon = self.DxH.dot(Mlon.dot(self.Dxu))
    Mbas = scipy.sparse.diags([basalCoeff],[0],shape=(self.Nx,self.Nx), format='csr')
    M = Mlon + Mbas

    rhs = -driving_k
    rhs[0] = 0.
    rhs[-1] = calving_ux # 0.5*self.delta*Hk[-1]**2
    
    M = M.asformat('lil') # seems to be happier about inserting entries in lil format
    M[0,0] = 1.
    M[0,1] = 0
    
#    # ux(xg) = (delta*H(xg)/(8*epsilon))**n
#    for index in range(len(glIndices)):
#      glIndex = glIndices[index]
#      alpha = lambda_g[index]
#      M[glIndex,glIndex-1] = (1-alpha)*longiCoeff[glIndex]*self.Dxu[glIndex,glIndex-1]
#      M[glIndex,glIndex] = (1-alpha)*longiCoeff[glIndex]*self.Dxu[glIndex,glIndex] \
#         + alpha*longiCoeff[glIndex+1]*self.Dxu[glIndex+1,glIndex]
#      M[glIndex,glIndex+1] = alpha*longiCoeff[glIndex+1]*self.Dxu[glIndex+1,glIndex+1]
#      Hxg = (1-alpha)*self.Hf[glIndex] + alpha*self.Hf[glIndex+1]
#      rhs[glIndex] = 0.5*self.delta*Hxg**2

    #M[-1,-2] = longiCoeff[-1]*self.Dxu[-1,-2]
    #M[-1,-1] = longiCoeff[-1]*self.Dxu[-1,-1]
    M[-1,-2] = self.Dxu[-1,-2]
    M[-1,-1] = self.Dxu[-1,-1]
    M = M.asformat('csr')
    
    norm = numpy.amax(numpy.abs(driving_k))
    
    resTau_k /= norm
    #print 'tau res k:', numpy.amax(numpy.abs(resTau_k))
       
 
    resTau_k2 = (M.dot(uk) - rhs)/norm
    
    #diff = numpy.amax(numpy.abs(resTau_k-resTau_k2))
    #print 'tau res k diff:', diff
    #print 'tau M norm:', norm
   
    ukp1 = scipy.sparse.linalg.spsolve(M,rhs)
    
#    uxkp1 = self.Dxu.dot(ukp1)
#    for index in range(len(glIndices)):
#      glIndex = glIndices[index]
#      alpha = lambda_g[index]
#      ux_xg = (1-alpha)*longiCoeff[glIndex]*uxkp1[glIndex] + alpha*longiCoeff[glIndex+1]*uxkp1[glIndex+1]
#      print 'internal bc error:', ux_xg-rhs[glIndex], alpha

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
#      fig = plt.figure(5)
#      if(len(fig.axes) > 0):
#        fig.axes[0].cla()
#      indices = numpy.arange(glIndices[0]-10,glIndices[0]+11)
#      #print xu[glIndices[0]]
#      plt.plot(xu[indices], resTau_k[indices]*norm, 'k',
#               xu[indices], resTau_k2[indices]*norm, 'k--')
#      fig = plt.figure(5)
#      if(len(fig.axes) > 0):
#        fig.axes[0].cla()
#      indices = numpy.arange(glIndices[0]-10,glIndices[0]+11)
#      #print xu[glIndices[0]]
#      plt.plot(xu[indices], deltaLongi[indices], 'b',
#               xu[indices], deltaBasal[indices], 'r',
#               xu[indices], deltaDriving[indices], 'g')
      if(self.plotContinuous):
        plt.draw()

    
    return (ukp1,numpy.max(numpy.abs(resTau_k2)))
    
  def solveCont(self, Hk, ukp1, HPrev, deltaT, firstIter):
    a = self.a # non-dimensional accumulation
    xH = self.xH
    #xu = self.xu
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

    #resCont_k = float(self.transient)*(Hk-HPrev)/deltaT + self.timeCentering*self.Dxu.dot(Hu) \
    #  + (1-self.timeCentering)*self.prevHux - a
    #resCont_k[0] = Hxk[0]
    
    M = self.timeCentering*self.Dxu.dot(Mup) + float(self.transient)*Mdt
    
    M = M.asformat('lil') # seems to be happier about inserting entries in lil format
    M[0,0] = self.DxH[0,0]
    M[0,1] = self.DxH[0,1]
    M = M.asformat('csr')
    
    #norm = numpy.amax(numpy.abs(HPrev/deltaT))
    #resCont_k /= norm

    rhs = a + self.transient*HPrev/deltaT - (1-self.timeCentering)*self.prevHux
    rhs[0] = 0.
    
    #resCont_k2 = (M.dot(Hk) - rhs)/norm
    
    #diff = numpy.amax(numpy.abs(resCont_k-resCont_k2))
    #print 'cont diff:', diff
    
    Hkp1 = scipy.sparse.linalg.spsolve(M,rhs)
      
    floatingMaskUk = self.floatingMaskU
    floatingMaskHk = self.floatingMaskH
    glIndicesk = self.glIndices
    lambda_gk = self.lambda_g

    (self.groundedMaskH,self.floatingMaskH,self.groundedMaskU,self.floatingMaskU,
     self.glIndices,self.lambda_g) \
      = self.getFlotationMasks(Hkp1)

    #resCont_kp1 = (M.dot(Hkp1) - rhs)/norm
    #print 'cont solver res:', numpy.amax(numpy.abs(resCont_kp1))
    #self.noiseFloor = numpy.maximum(self.noiseFloor, numpy.amax(numpy.abs(resCont_kp1)))
    
    Hu = Hkp1*ukp1
    Hu2 = Hkp1[1:]*ukp1[0:-1]
    mask = ukp1[0:-1] < 0.
    Hu[0:-1][mask] = Hu2[mask]  #never happens in practice: u >= 0 everywhere

    resCont_kp1 = float(self.transient)*(Hkp1-HPrev)/deltaT + self.Dxu.dot(Hu) - a
    resCont_kp1[0] = Hxk[0]
  
    if(self.plot):
      (sk,sxk) = self.computeSx(Hk,floatingMaskHk,floatingMaskUk,glIndicesk,lambda_gk)
      (skp1,sxkp1) = self.computeSx(Hkp1,self.floatingMaskH,
          self.floatingMaskU,self.glIndices,self.lambda_g)
      fig = plt.figure(3)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xH,sk, 'b', xH, skp1, 'r', xH, self.sPrev, 'm', xH, -self.b, 'k', 
               xH, self.Hf-self.b, 'g',
               xH[floatingMaskHk],sk[floatingMaskHk]-Hk[floatingMaskHk], 'b',
               xH[self.floatingMaskH],skp1[self.floatingMaskH]-Hkp1[self.floatingMaskH], 'r')
      fig = plt.figure(4)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xH,Hk, 'b', xH, Hkp1, 'r', xH, HPrev, 'm')
      if(self.plotContinuous):
        plt.draw()
    
    return Hkp1
    
  def step(self, HGuess, uGuess, HPrev, deltaT):
    if(self.plotContinuous):
      plt.ion()
  

    Hk = HGuess
    uk = uGuess
    (self.groundedMaskH,self.floatingMaskH,self.groundedMaskU,self.floatingMaskU,
     self.glIndices,self.lambda_g) \
      = self.getFlotationMasks(Hk)

    (groundedMaskHPrev,floatingMaskHPrev,groundedMaskUPrev,
     floatingMaskUPrev,
     self.glIndicesPrev,self.lambda_gPrev) \
      = self.getFlotationMasks(HPrev)
    (self.sPrev,sxPrev) = self.computeSx(HPrev,floatingMaskHPrev,
        floatingMaskUPrev,self.glIndicesPrev,self.lambda_gPrev)
   
    for iterIndex in range(self.maxPicardIter):
      (groundedMaskH,floatingMaskH,groundedMaskU,floatingMaskU,glIndices,lambda_g) \
        = self.getFlotationMasks(Hk)

      for inner in range(self.maxPicardInner):
        (ukp1,resTau) = self.iteratePicardTau(Hk, uk)
        uk = ukp1

      Hkp1 = self.solveCont(Hk, ukp1, HPrev, deltaT, (iterIndex == 0))
        
      tolerance = numpy.maximum(self.tolerance,2*self.noiseFloor)
      #print 'resTau, tol:', resTau, tolerance
      if(self.plot and not self.plotContinuous):
        plt.ioff()
        plt.show()  
        
      #print resTau, tolerance
      #test residuals for convergence...
      if((resTau < tolerance) and iterIndex > 0):
        print 'Picard converged in %i iterations.'%iterIndex
        print 'resTau, tol:', resTau, tolerance
        if(self.plot):
          plt.ioff()
          plt.show()
        return (Hkp1,ukp1)

      Hk = Hkp1
      #uk = ukp1
      
    print 'Warning:Picard iteration did not converge in', self.maxPicardIter, 'iterations. resTau:', resTau, 'tol:', tolerance
    #raise ValueError('Picard did not converge in %i iterations.'%self.maxPicardIter)
    return (Hkp1,ukp1)
    
  def computeHGuess(self,xg):
    xH = self.xH
    xu = self.xu
    n = self.n
    a = self.a
    
    glIndex = numpy.argmin(numpy.abs(xH-xg))
    xg = xH[glIndex]
    #shelfMask = xH > xg
    
    
    Hxg = self.computeB(xg)/(1.0-self.delta)
    if(Hxg <= 0.):
      raise ValueError("Hfself(xg) <= 0. Cannot produce HGuess")
        
    b = self.computeB(xH)
    
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
      
#    plt.figure(1)
#    plt.plot(xH,H, 'b', xH, HShelf, 'r')        
#    plt.figure(2)
#    plt.plot(xu, u, 'b', xu, uShelf, 'r')
#    plt.show()    
    
  
    return (H,u)
