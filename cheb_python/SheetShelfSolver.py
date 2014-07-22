#!/usr/bin/python
import numpy
import scipy.optimize
import Chebyshev
import numpy.linalg as linalg
from pprint import pprint
import initH
import os

def findSteadyStateValue(value, name, solution):
  setattr(solution, name, value)
  return solution.findSteadyState()
  
def findSteadyStateLog(logValue, name, solution):
  setattr(solution, name, 10**logValue)
  return solution.findSteadyState()
  
class SheetShelfSolver:
  def __init__(self, options):
    self.options = options

    # put all the variables that can be optimized into self
    self.A = options.A
    self.C = options.C
    self.lambda_0 = options.lambda_0
    self.linearSlope = options.linearSlope
    self.W = options.W
    self.p = options.p
    self.xg = options.xg
    
    self.Nx = options.Nx
    if(options.poly):
      self.computeB = self.computeBPoly
    else:
      self.computeB = self.computeBLinear

    self.cheb = Chebyshev.Chebyshev(self.Nx)
    
    self.delta = 1.0 - options.rho_i/options.rho_w

    sPerY = 365.25*24.*3600. # number of seconds per year
    self.aBar = .3/sPerY # m.s-1 accumulation rate
    self.HBar = 1000. # m hight scaling factor
    self.xBar = 1000000. # m domain scaling factor
    self.uBar = self.aBar*self.xBar/self.HBar # m/s ice velocity scaling factor
    self.WBar = 10000. # m width scaling factor
    
    self.updateNondimConstants()
    
    if(os.path.exists(options.inFile)):
      self.readResults(options.inFile)
      

    if(options.plot):
      import matplotlib.pyplot as plt

    if(options.plot and options.plotContinuous):
      plt.ion()

    print "self:"
    pprint(self.__dict__)
    
  def updateNondimConstants(self):
    self.NBar = options.rho_i*options.g*self.HBar
    self.taudBar = self.NBar*self.HBar/self.xBar
    self.taubBar = self.C*self.uBar**(1./options.n)
    self.Kappa = (options.m_0*self.uBar)/(self.lambda_0*options.Ab)/self.NBar**options.n
    self.gamma = self.taubBar/self.taudBar

    self.taulBar = (self.A**(-1./options.n))*(self.uBar/self.xBar)**(1./options.n)*self.HBar/self.xBar
    self.epsilon = self.taulBar/(2*self.taudBar)

    self.tauWbar = self.HBar/self.WBar*(5*self.uBar/(2*self.A*self.WBar))**(1./options.n)
    self.omega = self.tauWbar/self.taudBar
    

  def findSteadyStateBrentQ(self):
    
    options = self.options
    paramToOptimize = options.paramToOptimize
    if((paramToOptimize == 'A')
       or (paramToOptimize == 'C')
       or (paramToOptimize == 'lambda_0')
       or (paramToOptimize == 'linearSlope')
       or (paramToOptimize == 'W')):
      function = findSteadyStateLog
      lower = numpy.log10(options.minValue)
      upper = numpy.log10(options.maxValue)
    elif((paramToOptimize == 'xg')
       or (paramToOptimize == 'p')):
      function = findSteadyStateValue
      lower = options.minValue
      upper = options.maxValue
    else:
      raise ValueError("paramToOptimize not recognized")
      
  
    self.brentQIter = 0
    self.initWithPrev = True
    retry = False
    #try:
    value = numpy.array(scipy.optimize.brentq(function, lower, upper, 
                                             xtol=options.tolerance, args=(paramToOptimize,self)))
#    except ValueError as e:
#      print "Brentq failed:", e
#      if(self.brentQIter > 0):
#        retry = True
#      else:
#        return False

    if(retry or (numpy.abs(self.resBC) > options.tolerance)):
      print "Brentq failed: resBC is larger than the tolerance."
      print "Trying again without initializing from previous H"
      self.brentQIter = 0
      self.initWithPrev = False
      try:
        value = numpy.array(scipy.optimize.brentq(function, lower, upper, 
                                               xtol=options.tolerance, args=(paramToOptimize,self)))
      except ValueError as e:
        print "Brentq failed:", e
        return False

      if(numpy.abs(self.resBC) > options.tolerance):
        print "Brenq failed: resBC is larger than the tolerance."
        return False
    
    # make sure the the final solution is with the argument that minimizes the function
    function(value,paramToOptimize,self)
    
    return True

  def findSteadyState(self):
    options = self.options
    if(self.computeB(self.xg) <= 0.):
      self.resBC = 1e8
      raise ValueError("b <= 0 at grounding line. Cannot proceed with Brentq")
      
    print "BrentQ iteration: ", self.brentQIter+1
    print options.paramToOptimize, ":", getattr(self,options.paramToOptimize)
    
    self.initializeH()
  
    fracTolerance = 1e-5
    prevResBC = 0.0
    kMax = 200
    self.HScale = 1.0
    for k in range(kMax):
      print "k = ", k
      self.iterateUH()
      if(numpy.any(self.H  <= 0.)):
        self.resBC = 1e8
        raise ValueError("H <= 0 somewhere. Cannot proceed with Brentq")
  
      deltaResBC = numpy.abs(self.resBC-prevResBC)
      prevResBC = self.resBC
      if(self.resBC == 0.0):
        frac = 0.0
      else:
        frac = deltaResBC/numpy.abs(self.resBC)
      print "resBC:", self.resBC, deltaResBC, frac, fracTolerance
      print "resStress:", self.resStress, max(2*self.noiseFloor,options.tolerance)
      if((frac < fracTolerance) and (self.HScale == 1.0) and (numpy.abs(self.resBC) > options.tolerance)):
        # resBC isn't very close to zero, and it seems to be steady so there's no point if being too exact.
        break
      if(self.resStress < max(2*self.noiseFloor,options.tolerance)):
        break
    if(k == kMax-1):
      print "Warning: max iterations reached!"
    self.brentQIter += 1
    print "Picard iterations: ", k+1
    print "BrentQ iteration: ", self.brentQIter
    print "resBC", self.resBC
    
    if(self.initWithPrev):
      if(self.brentQIter > 1):
        self.Hprev2 = self.Hprev1
        self.xgprev2 = self.xgprev1
      self.Hprev1 = self.H
      self.xgprev1 = self.xg
    
    return self.resBC

  def initializeH(self):
    initialized = True
    (b_xg,bx_xg) = self.computeB(self.xg)
    Hf_xg = b_xg/(1.0-self.delta)
    if(self.initWithPrev and hasattr(self,'Hprev1')):
      if(self.brentQIter > 1):
        d1 = numpy.abs(self.xgprev1-self.xg)
        d2 = numpy.abs(self.xgprev2-self.xg)
        if(d1 < d2):
          H = self.Hprev1
        else:
          H = self.Hprev2
      else:
        H = self.Hprev1
    elif(hasattr(self,'initH')):
      H = self.initH
    else:
      initialized = False
      
    if(initialized):
      # we just need to scale the existing H to have H=Hf at the grounding line
      H0 = H[0]
      Hf_old = H[self.Nx-1]
      a = (H0-Hf_xg)/(H0-Hf_old)
      b = H0*(1. - a)
      H = a*H+b
      self.initH = H
    else:
      options = self.options
#      if(options.paramToOptimize == 'xg'):
#        (HGuessSheet,xgGuess) = initH.initHBounded(self,options.minValue,options.maxValue)
#        self.xg = xgGuess
#      else:
      (HGuessSheet,xgGuess) = initH.initH(self,self.xg)
      if(xgGuess > options.xc):
        print "xgGuess=%f > xc=%f: cannot proceed!"%(xgGuess,options.xc)
        exit(1)
      print "initial xg:", xgGuess
      H = numpy.zeros((2*self.Nx))
      H[0:self.Nx] = HGuessSheet
      uxg = options.a*xgGuess/HGuessSheet[-1]
      xShelf = xgGuess + (options.xc-xgGuess)*self.cheb.x
      const = (self.delta*options.a/(8.*self.epsilon))**options.n
      uShelf = (const*(xShelf**(options.n+1)-xgGuess**(options.n+1)) + uxg**(options.n+1))**(1./(options.n + 1))
      H[self.Nx:2*self.Nx] = options.a*xShelf/uShelf
      
      self.initH = H
      self.initXg = xgGuess
      self.H = H        

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

  def readResults(self):
    fileName = self.options.inFile
    print "reading from: ", fileName
    try:
      filePointer = open(fileName,'rb')
      Nx = numpy.fromfile(filePointer, dtype=int, count=1)[0]
      if(Nx != self.Nx):
        print "Bad Nx in file."
        exit(1)
  
      self.initXg = numpy.fromfile(filePointer, dtype=float, count=1)[0]
      self.initH = numpy.fromfile(filePointer, dtype=float, count=Nx)
      filePointer.close()
    except IOError:
      print 'Could not read %s. Exiting.'%fileName
      exit(40)
      
  def writeResults(self):
    fileName = "%s/%s"%(self.options.folder,self.options.outFile)
    print "writing to: ", fileName
    filePointer = open(fileName,'wb')
    Nx = numpy.array(self.Nx,int)
    Nx.tofile(filePointer)
    xg = numpy.array(self.xg)
    xg.tofile(filePointer)
    self.H.tofile(filePointer)
    self.u.tofile(filePointer)
    self.x.tofile(filePointer)
    value = numpy.array(self.A)
    value.tofile(filePointer)
    value = numpy.array(self.C)
    value.tofile(filePointer)
    value = numpy.array(self.lambda_0)
    value.tofile(filePointer)
    value = numpy.array(self.linearSlope)
    value.tofile(filePointer)
    value = numpy.array(self.W)
    value.tofile(filePointer)
    value = numpy.array(self.p)
    value.tofile(filePointer)
    filePointer.close()
 
  def absS(self,x):
    return numpy.sqrt(x**2 + self.options.eps_s**2)
    
  def absSx(self,x):
    return x/self.absS(x)
    
  def maxS(self,x):
    return 0.5*(self.absS(x)+x)
    
  def maxSx(self,x):
    return 0.5*(self.absSx(x) + 1.0)
    
  def derivX(self,field):
    options = self.options
    xg = self.xg
    Nx = self.Nx
    xc = options.xc
    DxSheet = self.cheb.Dx/xg
    DxShelf = self.cheb.Dx/(xc-xg)
    field_x = numpy.zeros(field.shape)
    field_x[0:Nx] = numpy.dot(DxSheet,field[0:Nx])
    field_x[Nx:2*Nx] = numpy.dot(DxShelf,field[Nx:2*Nx])
    return field_x
    
  def updateX(self):
    options = self.options
    xg = self.xg
    Nx = self.Nx
    xc = options.xc
    xSheet = xg*self.cheb.x
    xShelf = xg + (xc-xg)*self.cheb.x
    x = numpy.zeros((2*Nx))
    x[0:Nx] = xSheet
    x[Nx:2*Nx] = xShelf
    self.x = x 

  def iterateUH(self):
    
    self.updateNondimConstants()
    self.updateX()
 
    options = self.options
    xg = self.xg
    Nx = self.Nx
    xc = options.xc
    DxSheet = self.cheb.Dx/xg
    DxShelf = self.cheb.Dx/(xc-xg)
    Dx = numpy.zeros((2*Nx,2*Nx))
    Dx[0:Nx,0:Nx] = DxSheet
    Dx[Nx:2*Nx,Nx:2*Nx] = DxShelf
  
    p = self.p
    Kappa = self.Kappa
    gamma = self.gamma
    epsilon = self.epsilon
    omega = self.omega
    delta = self.delta
    n = options.n
    
    a = options.a
    W = self.W

    if(xg >= xc):
      print "Error: xg=%f >= xc=%f. Exiting..."%(xg,xc)
      exit(1)
      
    H = self.H
    x = self.x
    
    Hx = self.derivX(H)

    (b,bx) = self.computeB(x)
    
    Hf = self.maxS(b[0:Nx])/(1-delta)
    Np = H[0:Nx]*(self.maxS(1.0 - Hf/H[0:Nx]))**p

    s = numpy.zeros(x.shape)
    s[0:Nx] = H[0:Nx]-b[0:Nx]
    s[Nx:2*Nx] = delta*H[Nx:2*Nx]

    # driving stress = drivingCoeff1 Hx + drivingCoeff2
    drivingCoeff1 = numpy.zeros(x.shape)
    drivingCoeff2 = numpy.zeros(x.shape)
    drivingCoeff1[0:Nx] = -H[0:Nx]
    drivingCoeff2[0:Nx] = -H[0:Nx]*bx[0:Nx]
    drivingCoeff1[Nx:2*Nx] = -delta*H[Nx:2*Nx]
    
    u = a*x/H
    ux = (a - u*Hx)/H
    Hxx = self.derivX(Hx)
    uxx = (-2*Hx*ux - Hxx*u)/H
    abs_ux = self.absS(ux)
    nu = 4.*epsilon*(abs_ux)**(1./n-1.)
    abs_u = self.absS(u)
    basal = numpy.zeros(u.shape)
    if(options.useSchoofBasal):
      gammaPrime = self.C*self.uBar**(1./n)/(options.rho_i*options.g*self.HBar**2/self.xBar)
      basal[0:Nx] = -gammaPrime*abs_u[0:Nx]**(1./n)*u[0:Nx]/abs_u[0:Nx]
    else:
      basal[0:Nx] = -gamma*Np*(abs_u[0:Nx]/(Kappa*abs_u[0:Nx] + Np**n))**(1./n)*u[0:Nx]/abs_u[0:Nx]
    
    nux = nu*(1./n - 1)*uxx/abs_ux
    
    longiCoeff1 = nux*a
    longiCoeff2 = -(nux*u + nu*ux)
    longiCoeff3 = -nu*u
    longi = longiCoeff1 + longiCoeff2*Hx + longiCoeff3*Hxx
    
    driving = drivingCoeff1*Hx + drivingCoeff2
    
    
    #tau_w = -omega*H*W^(-1-1/n)*|u|^(1/n-1) u
    if(options.useChannelWidth):
      lateralCoeff = -omega*W**(-1.-1./n)*abs_u**(1/n-1)*u
    else:
      lateralCoeff = numpy.zeros(x.shape)
      
    lateral = lateralCoeff*H
    
    residual = longi + basal + driving + lateral
     
    # tau_l = H nu ux
    # ux = (a - u Hx)/H
    # tau_l = nu (a - u Hx)
    # tau_lx = nux a - (nu u Hx)_x
    #        = nux a - (nux u + nu ux) Hx - nu u Hxx
    
    # tau_d = -H sx
    
    Dxx = numpy.dot(Dx,Dx)
    
    M = numpy.dot(numpy.diag(longiCoeff2 + drivingCoeff1),Dx) \
      + numpy.dot(numpy.diag(longiCoeff3),Dxx) \
      + numpy.diag(lateralCoeff)
    rhs = -longiCoeff1 - basal - drivingCoeff2
    
    # Hx = bx at the ice divide
    M[0,:] = 0
    M[0,0:Nx] = DxSheet[0,:]
    rhs[0] =  bx[0]
    
    # H = Hf at the grounding line
    M[Nx-1,:] = 0.
    M[Nx-1,Nx-1] = 1.
    rhs[Nx-1] = Hf[Nx-1]
    M[Nx,:] = 0.
    M[Nx,Nx] = 1.
    rhs[Nx] = Hf[Nx-1]
    
    # Hx is continuous at the grounding line
    M[2*Nx-1,0:Nx] = DxSheet[Nx-1,:]
    M[2*Nx-1,Nx:2*Nx] = -DxShelf[0,:]
    rhs[2*Nx-1] = 0.
    
    
    newH = linalg.solve(M,rhs)
    
    solver_res = numpy.dot(M,newH) - rhs
    residual_check = numpy.dot(M,H) - rhs
    self.noiseFloor = numpy.maximum(numpy.max(numpy.abs(solver_res[1:Nx-1])),
                                    numpy.max(numpy.abs(solver_res[Nx+1:2*Nx-1])))
                                    
                                    
    minScale = 0.1
    self.HScale = 1.0
    if(numpy.any(newH < minScale*H)):
      print "Have to scale newH so H doesn't got through zero." 
      self.HScale = numpy.amin(newH/H)
      print self.HScale
      alpha = (1.0-minScale)/(1.0-self.HScale)
      print alpha
      newH = alpha*newH + (1-alpha)*H

    newU = a*x/newH  
    
    if(options.plot):
      import matplotlib.pyplot as plt
      newS = numpy.zeros(x.shape)
      newS[0:Nx] = newH[0:Nx]-b[0:Nx]
      newS[Nx:2*Nx] = delta*newH[Nx:2*Nx]
      if(options.plotContinuous):
        plt.ion()
      temp = nu[-1]*ux[Nx-1]/(0.5*self.delta) - b[Nx-1]
      fig = plt.figure(1)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(x,s, 'b', x, newS, 'r', x, -b, 'k', x[0:Nx], Hf-b[0:Nx], 'g',
               x, temp*numpy.ones(x.shape),'k--')
      fig = plt.figure(2)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(x,solver_res,'b', x, residual, 'r',
               x, residual_check, 'k--')
      fig = plt.figure(3)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(x,longi, 'b', x, basal, 'r', x, driving, 'g', x, residual, 'k',
               x, residual_check, 'k--')
      fig = plt.figure(4)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(x,u, 'b', x, newU, 'r', x, a*xg/Hf[Nx-1]*numpy.ones(x.shape,float),'--k')
      fig = plt.figure(5)
      if(len(fig.axes) > 0):
        fig.axes[0].cla()
      plt.plot(x,ux,'b',x,a/H,'r',x,-u*Hx/H,'g')
      if(options.plotContinuous):
        plt.draw()
        plt.pause(0.0001)
      else:
        plt.show()
      

    self.resBC = nu[Nx-1]*ux[Nx-1] - 0.5*delta*Hf[Nx-1]
    
    self.resStress = numpy.maximum(numpy.max(numpy.abs(residual_check[1:Nx-1])),
                                   numpy.max(numpy.abs(residual_check[Nx+1:2*Nx-1])))
    self.H = H
    self.u = newU
    
from optparse import OptionParser
    
parser = OptionParser()

parser.add_option("--minValue", type="float", default=1.0, dest="minValue")
parser.add_option("--maxValue", type="float", default=2.0, dest="maxValue")
parser.add_option("--paramToOptimize", type="string", default="xg", dest="paramToOptimize")

parser.add_option("--xc", type="float", default=2.112, dest="xc")
parser.add_option("--useChannelWidth", action="store_true", dest="useChannelWidth")
parser.add_option("--W", type="float", default=1.0, dest="W")

parser.add_option("--xg", type="float", default=1.0, dest="xg")
parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=2.1544e-25, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--rho_i", type="float", default=900.0, dest="rho_i")
parser.add_option("--rho_w", type="float", default=1000.0, dest="rho_w")
parser.add_option("--g", type="float", default=9.8, dest="g")
parser.add_option("--n", type="float", default=3.0, dest="n")
parser.add_option("--a", type="float", default=1.0, dest="a")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007
parser.add_option("--lambda_0", type="float", default=2.0, dest="lambda_0")
parser.add_option("--m_0", type="float", default=0.5, dest="m_0")
parser.add_option("--Ab", type="float", default=3.1688e-24, dest="Ab")
parser.add_option("--poly", action="store_true", dest="poly")

parser.add_option("--inFile", type="string", default="none", dest="inFile")
parser.add_option("--folder", type="string", default="results", dest="folder")
parser.add_option("--outFile", type="string", default="results.pyda", dest="outFile")

parser.add_option("--Nx", type="int", default=513, dest="Nx")

parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")
parser.add_option("--tolerance", type="float", default=1e-10, dest="tolerance")

parser.add_option("--plot", action="store_true", dest="plot")
parser.add_option("--plotContinuous", action="store_true", dest="plotContinuous")
parser.add_option("--useSchoofBasal", action="store_true", dest="useSchoofBasal")

options, args = parser.parse_args()

solver = SheetShelfSolver(options)


success = solver.findSteadyStateBrentQ()
if(not success):
  print "Did not find a solution, not writing to a file."
  exit()

solver.writeResults()