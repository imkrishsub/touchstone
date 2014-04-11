#!/usr/bin/python
import numpy
#import math
import scipy.optimize
import Chebyshev
import numpy.linalg as linalg
import matplotlib.pyplot as plt

def findSteadyState_A_xg(A,xg, solution):
  if(solution.computeB(xg,solution) <= 0.):
    solution.resBC = 1e8
    raise ValueError("b <= 0 at grounding line. Cannot proceed with Brentq")
    
  solution.xg = numpy.array(xg)
  solution.A = numpy.array(A)
  print "xg = ", xg
  print " p = ", solution.p
  print " A = ", A
  print "BrentQ iteration: ", solution.brentQIter+1
  
  x = solution.cheb.x*xg
  if(solution.plotContinuous):
    plt.ion()

  (solution.H, solution.Hx) = solution.computeHGuess(x,solution)

  prevResBC = 0.0
  resBCTolerance = 1e-3
  kMax = 200
  for k in range(kMax):
    print "k = ", k
    solution.iterateUH()
    if(numpy.any(solution.H  <= 0.)):
      solution.resBC = 1e8
      raise ValueError("H <= 0 somewhere. Cannot proceed with Brentq")

    deltaResBC = numpy.abs(solution.resBC-prevResBC)
    prevResBC = solution.resBC
    if(solution.resBC == 0.0):
      frac = 0.0
    else:
      frac = deltaResBC/numpy.abs(solution.resBC)
    print "resBC:", solution.resBC, deltaResBC, frac, resBCTolerance
    print "resStress:", solution.resStress, max(2*solution.noiseFloor,solution.tolerance)
    if((frac < resBCTolerance) and (numpy.abs(solution.resBC) > solution.tolerance)):
      # resBC isn't very close to zero, and it seems to be steady so there's no point if being too exact.
      break
    if(solution.resStress < max(2*solution.noiseFloor,solution.tolerance)):
      break
  if(k == kMax-1):
    print "Warning: max iterations reached!"
  solution.brentQIter += 1
  print "A =  ", solution.A
  print "p =  ", solution.p
  print "xg = ", solution.xg
  print "Picard iterations: ", k+1
  print "BrentQ iteration: ", solution.brentQIter
  print "resBC", solution.resBC
  
  if(solution.initWithPrev):
    if(solution.brentQIter > 1):
      solution.Hprev2 = solution.Hprev1
      solution.Hxprev2 = solution.Hxprev1
      solution.xgprev2 = solution.xgprev1
    solution.Hprev1 = solution.H
    solution.Hxprev1 = solution.Hx
    solution.xgprev1 = solution.xg
  
  return solution.resBC

def findSteadyState_xg(xg, solution):
  return findSteadyState_A_xg(solution.A,xg, solution)
  
def findSteadyState_A(logA, solution):
  solution.logA = logA
  return findSteadyState_A_xg(10**logA,solution.xg, solution)
  
class SteadyStateSolution:
  def __init__(self, Nx):
    self.Nx = Nx
    self.cheb = Chebyshev.Chebyshev(Nx)
    self.intX_H = self.cheb.intX.copy()
    for index in range(Nx):
      self.intX_H[index,:] -= self.intX_H[-1,:]
    self.rho_i = 900.0 # kg/m^3 ice density
    self.rho_w = 1000.0 # kg/m^3 water density
    self.delta = 1.0 - self.rho_i/self.rho_w
    self.C = numpy.array(7.624e6)
    self.a = numpy.array(1.0)
    self.plotContinuous = True
    self.plot = True
    self.useLongi = True
    self.useSchoofBasal = False
    #self.Ab = 4.227e-25
    self.Ab = 3.1688e-24 #Pa^-3 s^-1
    
    sPerY = 365.25*24.*3600. # number of seconds per year
    self.g = 9.8 # m/s^2 gravity acceleration
    self.aBar = .3/sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    self.n = 3. # Glen's flow parameter

    self.lambda_0 = 2   # wavelength of bedrock bump (m)
    self.m_0  = .5  # maximum bed obstacle slope (no unit)
    
    
    
  def findSteadyStateBrentQxg(self, p, A, xgMin, xgMax, eps_s, tolerance, computeB, computeHGuess):
    self.A = numpy.array(A)
    self.p = numpy.array(p)
    self.eps_s = numpy.array(eps_s)
    self.tolerance = tolerance
    self.computeB = computeB
    self.computeHGuess = computeHGuess
    self.alpha = 0.0
  
    (success, xgMin, xgMax) = self.findBrentqLimits(findSteadyState_xg, xgMin, xgMax)
    
    if(not success):
      print "Brentq bounds not found."
      return False
      
        
    self.brentQIter = 2
    self.initWithPrev = True
    try:
      xg = numpy.array(scipy.optimize.brentq(findSteadyState_xg, xgMin, xgMax, xtol=tolerance, args=(self,)))
    except ValueError as e:
      print "Brentq failed:", e
      print "Trying again without initializing from previous Hx"
      self.brentQIter = 0
      self.initWithPrev = False
      try:
        xg = numpy.array(scipy.optimize.brentq(findSteadyState_xg, xgMin, xgMax, xtol=tolerance, args=(self,)))
      except ValueError as e:
        print "Brentq failed:", e
        return False

    if(numpy.abs(self.resBC) > self.tolerance):
      print "Brenq failed: resBC is larger than the tolerance."
      return False
      
    if(self.xg != xg):
      self.xg = xg
      findSteadyState_xg(self.xg, self)
    return True

  def findSteadyStateBrentQA(self, p, xg, AMin, AMax, eps_s, tolerance, computeB, computeHGuess):
    #self.A = numpy.array(A)
    self.xg = numpy.array(xg)
    self.p = numpy.array(p)
    self.eps_s = numpy.array(eps_s)
    self.tolerance = tolerance
    self.computeB = computeB
    self.computeHGuess = computeHGuess
    self.alpha = 0.0

    (success, logAMin, logAMax) = self.findBrentqLimits(findSteadyState_A, numpy.log10(AMin), numpy.log10(AMax))
    
    if(not success):
      print "Brentq bounds not found."
      return False

    self.brentQIter = 2

    self.initWithPrev = True
    try:
      logA = numpy.array(scipy.optimize.brentq(findSteadyState_A, logAMin, logAMax, xtol=tolerance, args=(self,)))
    except ValueError as e:
      print "Brentq failed:", e
      return False
    
    if(numpy.abs(self.resBC) > self.tolerance):
      print "Brenq failed: resBC is larger than the tolerance."
      return False
    if(self.logA != logA):
      findSteadyState_A(logA, self)
    return True

  def findSteadyState(self, p, A, xg, eps_s, tolerance, computeB, computeHGuess):
    self.A = numpy.array(A)
    self.p = numpy.array(p)
    self.eps_s = numpy.array(eps_s)
    self.computeB = computeB
    self.computeHGuess = computeHGuess
    self.brentQIter = 0
    self.tolerance = tolerance
    
    self.alpha = 1.0 
    findSteadyState_xg(xg,self)
    
    
  def findBrentqLimits(self, fun, minValueOrig, maxValueOrig):
  
    self.initWithPrev = False

    minValueOuter = minValueOrig
    maxValueOuter = maxValueOrig
    
    minValue = minValueOuter
    maxValue = maxValueOuter
    success = False


    for iteration in range(10):
      print "iteration", iteration, minValue, maxValue, minValueOuter, maxValueOuter
      try:
        self.brentQIter = 0
        self.plot = True
        fa = fun(minValue, self)
        print "fa:", fa
        self.Hprev2 = self.H
        self.Hxprev2 = self.Hx
        self.xgprev2 = self.xg
        minValid = True
      except ValueError as e:
        print "Error:", e
        minValid = False
        # increase minValue
    
  
      try:
        self.brentQIter = 0
        self.plot = True
        fb = fun(maxValue, self)
        print "fb:", fb
        self.Hprev1 = self.H
        self.Hxprev1 = self.Hx
        self.xgprev1 = self.xg
        maxValid = True
      except ValueError as e:
        print "Error:", e
        # decrease maxValue
        maxValid = False
        
      alphaMin = 0.5
      alphaMax = 0.5

      if(not minValid):
        minValueOuter = minValue 
      if(not maxValid):
        maxValueOuter = maxValue
      
      if(not minValid and not maxValid):
        alphaMin = 0.25
        alphaMax = 0.25

      if(not minValid):
        minValue = (1.0-alphaMin)*minValueOuter+alphaMin*maxValueOuter
      if(not maxValid):
        maxValue = alphaMax*minValueOuter+(1.0-alphaMax)*maxValueOuter

      if(minValid and maxValid):
        success = (numpy.sign(fa) != numpy.sign(fb))
        if(success):
          break
        
        if((minValue == minValueOuter) and (maxValue == maxValueOuter)):
          # can't expand so we've failed to find valid bounds
          print "Could not find valid Brentq bounds: f(a) has the same sign as f(b)"
          break
        
        # expand 
        minValue = 0.5*(minValue + minValueOuter)
        maxValue = 0.5*(minValue + minValueOuter)
        print "min, max:", minValue, maxValue
  
    if(not success and iteration == 9):
      print "Did not find valid bounds in 10 iterations."
      
    if(success):
      print "Success! Ready to proceed with Brentq with bounds", minValue, maxValue
    return (success, minValue, maxValue)


  
  def absS(self,x):
    return numpy.sqrt(x**2 + self.eps_s**2)
    
  def absSx(self,x):
    return x/self.absS(x)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)
    
  def maxSx(self,x):
    return 0.5*(self.absSx(x) + 1.0)

  def iterateUH(self):
    xg = self.xg
    Dx = self.cheb.Dx/xg
    #intX_H = self.intX_H*xg
    xk = xg*self.cheb.x
    p = self.p
    A = self.A
#    fig = plt.figure(1)
#    fig.add_subplot(111)
#    plt.plot(x,numpy.dot(Dx,x**2))
#    plt.show()
#    crash
    
    
    n = self.n # Glen's flow parameter
    m = 1./n  

    ####################################################################
    
                          # Basal friction constants #
    
                          
    # values taken from Pimentel & al (2010)                      
    
    #Lambda_0 = self.lambda_0*A/self.m_0 
    Lambda_0 = self.lambda_0*self.Ab/self.m_0 
        
    C_tilda = self.C*Lambda_0**m    
        
    Kappa = (Lambda_0*(self.rho_i*self.g*self.HBar)**n)/self.uBar # non-dimensionalization coeff in 2nd term in tau_b
    
    ####################################################################

    epsilon = ((A**(-1./n))*(self.uBar/self.xBar)**(1./n))/(2*self.rho_i*self.g*self.HBar)
    #if(self.useLongi):
    #    epsilon = ((A**(-1./n))*(self.uBar/self.xBar)**(1./n))/(2*self.rho_i*self.g*self.HBar)
    #else:
    #    epsilon = 0
    
    gamma = C_tilda*self.xBar/self.HBar
    
#    print gamma/Kappa**(1/n)
#    print 1/Kappa
#    print epsilon
#    crash
    
    a = self.a # non-dimensional accumulation
    
    bk = self.computeB(xk,self)
    bxk = numpy.dot(Dx,bk) 
    Hfk = self.maxS(bk)/(1-self.delta)

    Hxk = self.Hx
    #Hk = numpy.dot(intX_H,Hxk) + Hfk[-1] #Chebyshev.intXPhys(Hxk,xg) 
    Hk = Chebyshev.intXPhys(Hxk,xg) 
    Hk = Hk + (Hfk[-1] - Hk[-1])
    uk = a*xk/Hk
    uxk = (a - uk*Hxk)/Hk
    Hxxk = numpy.dot(Dx,Hxk)
    uxxk = (-2*Hxk*uxk - Hxxk*uk)/Hk
    abs_ux = self.absS(uxk)
    nuk = 4.*epsilon*(abs_ux)**(1./n-1.)
    fpk = (self.maxS(1.0 - Hfk/Hk))**p
    Npk = Hk*fpk
    abs_u = self.absS(uk)
    if(self.useSchoofBasal):
      gammaPrime = self.C*self.uBar**m/(self.rho_i*self.g*self.HBar**2/self.xBar)
      #basal = -gamma*(abs_u/Kappa)**(1./n)*uk/abs_u
      basal = -gammaPrime*abs_u**(1./n)*uk/abs_u
    else:
      basal = -gamma*Npk*(abs_u/(abs_u + Kappa*Npk**n))**(1./n)*uk/abs_u
    
    nuxk = nuk*(1./n - 1)*uxxk/abs_ux
    if(self.useLongi):
      longi = nuxk*a -(nuxk*uk + nuk*uxk)*Hxk - nuk*uk*Hxxk
    else:
      longi = numpy.zeros(xk.shape)
    #longi_f = Chebyshev.transformPhysToCheb(longi)
    

    #diff = bk-(1-self.delta)*Hk
    #diffx = bxk-(1-self.delta)*Hk
    # absSx = x/absS(x)
    # maxSx = 0.5*(absSx + 1)

    #botk = -bk + self.maxS(diff)
    #botxk = -bxk + self.maxSx(diff)*diffx
    #sk = botk + Hk
    #sxk = botxk + Hxk
    #sxCheck = numpy.dot(Dx,sk)
    
    #fig = plt.figure(1)
    #plt.plot(xk,sxk, 'b', xk, sxCheck, 'r')
    #plt.show()

    sk = Hk-bk
    sxk = Hxk-bxk
    
    driving = -Hk*sxk #-Hk*(Hxk - bxk)
    
    residual = longi + basal + driving
    
    
    #if(numpy.any(Hk < Hfk)):
    #  raise ValueError("H < Hf")
    #if(numpy.any(uxk < 0)):
    #  raise ValueError("ux < 0")
    # H = numpy.dot(intX,Hx)
    # H = H - H[-1] + Hf[-1]
    #   = numpy.dot(intX,Hx) - numpy.dot(intX[-1,:],Hx) + Hf[-1]
    #   = numpy.dot(intX - intX[-1,:],Hx) + Hf[-1]

    #M = -numpy.diag(nuxk*uk + nuk*uxk + Hk) - numpy.dot(numpy.diag(nuk*uk),Dx)
    #rhs = -nuxk*a - basal - Hk*bxk
    if(self.useLongi):
      M = -numpy.diag(nuxk*uk + nuk*uxk + Hk) - numpy.dot(numpy.diag(nuk*uk),Dx)
      rhs = -nuxk*a - basal - Hk*bxk #Hk*(sxk-Hxk)
    else:
      M = -numpy.diag(Hk)
      rhs = - basal - Hk*bxk #Hk*(sxk-Hxk)
    M[0,:] = 0.
    M[0,0] = 1.
    rhs[0] =  0.
    Hxkp1 = linalg.solve(M,rhs)
    Hkp1 = Chebyshev.intXPhys(Hxkp1,xg)
    Hkp1 = Hkp1 + (Hfk[-1] - Hkp1[-1])


#    mask = Hkp1 < Hfk
#    Hkp1[mask] = Hfk[mask]
#    Hfxk = numpy.dot(Dx,Hfk)
#    Hxkp1[mask] = Hfxk[mask]
    
    
    solver_res = numpy.dot(M,Hxkp1) - rhs
    residual_check = numpy.dot(M,Hxk) - rhs
    self.noiseFloor = numpy.max(numpy.abs(solver_res[1:-1]))

    minScale = 0.1
    if(numpy.any(Hkp1 < minScale*Hk)):
      print "Have to scale Hkp1 so H doesn't got through zero." 
      scale = numpy.amin(Hkp1/Hk)
      print scale
      alpha = (1.0-minScale)/(1.0-numpy.amin(scale))
      Hxkp1 = alpha*Hxkp1 + (1-alpha)*Hxk
      Hkp1 = alpha*Hkp1 + (1-alpha)*Hk

    ukp1 = a*xk/Hkp1
    uxkp1 = (a - ukp1*Hxkp1)/Hkp1

    if(numpy.any(uxkp1 < minScale*uxk)):
      print "Have to scale uxkp1 so ux doesn't got through zero." 
      scale = numpy.amin(uxkp1/uxk)
      print scale
      alpha = (1.0-minScale)/(1.0-numpy.amin(scale))
      print alpha
      Hxxkp1 = numpy.dot(Dx,Hxkp1)
      uxxkp1 = (-2*Hxkp1*uxkp1 - Hxxkp1*ukp1)/Hkp1
      ukp1 = alpha*ukp1 + (1-alpha)*uk
      uxkp1 = alpha*uxkp1 + (1-alpha)*uxk
      uxxkp1 = alpha*uxxkp1 + (1-alpha)*uxxk
      Hkp1[1:] = a*xk[1:]/ukp1[1:]
      Hkp1[0] = a/uxkp1[0]
      Hxkp1[1:] = (a - Hkp1[1:]*uxkp1[1:])/ukp1[1:]
      Hxkp1[0] = -0.5*(uxxkp1[0]*Hkp1[0])/uxkp1[0]

    
    #if(numpy.any(uxkp1 < 0.0)):
      #raise ValueError("Iteration Cannot proceed -- getting ux < 0")
      
#    diffk = Hk - xk*Hxk
#    diffkp1 = Hkp1 - xk*Hxkp1
#    
#    minDiff = 0.1
#    limiter = (1.0 - minDiff)*diffk/(diffk-diffkp1)
#    mask = diffkp1 >= minDiff*diffk
#    limiter[mask] = 1.0

    
    diff = bk-(1-self.delta)*Hkp1
    botkp1 = -bk + self.maxS(diff)
    skp1 = botkp1 + Hkp1

    if(self.plot):
      temp = nuk[-1]*uxk[-1]/(0.5*self.delta) - bk[-1]
      fig = plt.figure(1)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xk,sk, 'b', xk, skp1, 'r', xk, -bk, 'k', xk, Hfk-bk, 'g',
               #xk,botk, 'c', xk, botkp1, 'm',
               xk, temp*numpy.ones(xk.shape),'k--')
               #xk, xk*Hxkp1 - bk, 'm', xk, xk*Hxk - bk, 'c')
      fig = plt.figure(2)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xk,solver_res,'b', xk, residual, 'r',
               xk, residual_check, 'k--')
      fig = plt.figure(3)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xk,longi, 'b', xk, basal, 'r', xk, driving, 'g', xk, residual, 'k',
               xk, residual_check, 'k--', xk, -rhs, 'c',
               xk, nuxk*a, 'm', xk, Hk*bxk, 'y')
      fig = plt.figure(4)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(xk,uk, 'b', xk, ukp1, 'r', xk, a*xg/Hfk[-1]*numpy.ones(xk.shape,float),'--k')
      fig = plt.figure(5)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      #plt.semilogy(numpy.abs(longi_f))
      plt.plot(xk,uxk,'b',xk,a/Hk,'r',xk,-uk*Hxk/Hk,'g')
      fig = plt.figure(6)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      #plt.plot(xk,Hxk, 'b', xk, Hxkp1, 'r', xk, numpy.dot(Dx,Hk),'k--')
      plt.plot(xk,uxkp1/uxk,'b')
      if(self.plotContinuous):
        plt.draw()
      else:
        plt.show()
      
      
    if(self.alpha > 0.0):
      Hfx_xg = numpy.dot(Dx[-1,:],Hfk)
      print 'Hfx at xg: ', Hfx_xg
      xgkp1 = (nuk[-1]*uxk[-1]/(0.5*self.delta) - Hfk[-1])/Hfx_xg + xg
      self.xg = self.alpha*xgkp1 + (1.0-self.alpha)*xg 
      print "xg old, new, diff:", xg, xgkp1, xgkp1-xg
    
    self.resBC = nuk[-1]*uxk[-1] - 0.5*self.delta*Hfk[-1]
    
    self.resStress = numpy.max(numpy.abs(residual_check[1:-1]))
    


    self.Hx = Hxkp1
    self.H = Hkp1
    self.u = ukp1
    
    
# vim: ts=2 noet
