#!/usr/bin/python
import numpy
#import math
import matplotlib.pyplot as plt
import os.path
from SteadyStateSolution import SteadyStateSolution
from optparse import OptionParser


def extremaBPoly(solution):
  schoofx = 750000.
  #bx2 == 0 
  #    == -(-2184.8 + 2.0*1031.72*xs**2 - 3.0*151.72*xs**4)/solution.HBar 
  # xs**2 == (-b +- sqrt(b**2 - 4*a*c))/2a
  ra = -3.0*151.72
  rb = 2.0*1031.72
  rc = -2184.8
  radicand = rb**2 - 4.0*ra*rc
  
  xs2_r1 = (-rb + numpy.sqrt(radicand))/(2.0*ra)
  xs2_r2 = (-rb - numpy.sqrt(radicand))/(2.0*ra)
  
  x1 = schoofx/solution.xBar*numpy.sqrt(xs2_r1)
  x2 = schoofx/solution.xBar*numpy.sqrt(xs2_r2)
  return (x1, x2)


def computeBPoly(x,solution):
  schoofx = 750000.
  xs = x*solution.xBar/schoofx
  b = -(729 - 2184.8*xs**2 + 1031.72*xs**4 - 151.72*xs**6)/solution.HBar  
  return b

#def computeBxPoly(x):
#  schoofx = 750000.
#  HBar = 1000. # m hight scaling factor
#  xBar = 1000000. # m domain scaling factor
#  xs = x*xBar/schoofx
#  bx = -xBar/schoofx*(- 2.0*2184.8*xs + 4.0*1031.72*xs**3 - 6.0*151.72*xs**5)/HBar  
#  return bx

 
def computeHGuess(x,solution):
  H0 = 3.0
  Hf_xg = solution.computeB(x[-1],solution)/(1.0-solution.delta)
  #print Hf_xg
  if(Hf_xg <= 0.):
    raise ValueError("Hf(xg) <= 0. Cannot produce HGuess")
  if(solution.initWithPrev and hasattr(solution,'Hprev1')):
    if(solution.brentQIter > 1):
      d1 = numpy.abs(solution.xgprev1-solution.xg)
      d2 = numpy.abs(solution.xgprev2-solution.xg)
      if(d1 < d2):
        H = solution.Hprev1
        Hx = solution.Hxprev1*solution.xgprev1
      else:
        H = solution.Hprev2
        Hx = solution.Hxprev2*solution.xgprev2
    else:
      H = solution.Hprev1
      Hx = solution.Hxprev1*solution.xgprev1

    # we just need to scale the existing H to have H=Hf at the grounding line
    H0 = H[0]
    Hf_old = H[-1]
    a = (H0-Hf_xg)/(H0-Hf_old)
    b = H0*(1. - a)
    H = a*H+b
    Hx = a*Hx/x[-1]
  else:
    H = (Hf_xg - H0)*(x/x[-1])**2.0 + H0
    Hx = 2.0*(Hf_xg - H0)*(x/x[-1]**2)
  
  return (H,Hx)


parser = OptionParser()
parser.add_option("--region", type="int", default=1, dest="region", 
                  help="region number (1,2,3 -- retreated, advanced or unstable regions")

parser.add_option("--dense", action="store_true", dest="dense", 
                  help="densely sample A?")

parser.add_option("--useSchoofBasal", action="store_true", dest="useSchoofBasal", 
                  help="use the Schoof 2007 basal sliding law?")

parser.add_option("--noLongi", action="store_true", dest="noLongi", 
                  help="exclude longitudinal stress in stress balance?")

parser.add_option("--pIndex", type="int", default=0, dest="pIndex", 
                  help="index into p array (0, 0.25, 0.5, 0.75, 1)")

options, args = parser.parse_args()


Nx = 1025
#eps_s = 1e-3

if(options.dense):
  oldAs = [30,
    25,
    20,
    15,
    10,
    8.5,
    6.5,
    5,
    2.5,
    2,
    1.5,
    1,
    0.5,
    0.25]
  oldAs = numpy.array(oldAs)*1e-26 #Pa^(-3) s^(-1)
  
  density = 8
  As = numpy.zeros((density*len(oldAs)-(density-1)),float)
  As[0::density] = oldAs
  for offset in range(1,density):
    frac = offset/float(density)
    As[offset::density] = (1.0 - frac)*oldAs[0:-1]+frac*oldAs[1:] 
  #As = 10.0**numpy.linspace(-24.5,-26.6,41)
else:
  As = [30,
    25,
    20,
    15,
    10,
    5,
    2.5,
    2,
    1.5,
    1,
    0.5,
    0.25]
  As = numpy.array(As)*1e-26 #Pa^(-3) s^(-1)


eps_s = 1e-10
tolerance = 1e-10

if(options.useSchoofBasal):
  if(options.dense):
    root = './schoof_dense_Ab'
  else:
    root = './schoof_poly_Ab'
else:
  if(options.dense):
    root = './dense_Ab'
  else:
    root = './poly_Ab'


if (not os.path.exists(root)):
	os.mkdir(root)


ps = numpy.linspace(0.0,1.0,5)

#ps = ps[-1:]
#As = As[-2:]
#xgMax = 0.5*(xgMax+xgMin)

solution = SteadyStateSolution(Nx)

#import Chebyshev
#xMax = 1.5
#x = solution.cheb.x*xMax
#b = computeBPoly(x,solution)
#bx1 = computeBxPoly(x)
#bx2 = numpy.dot(solution.cheb.Dx,b)/xMax
#bx3 = Chebyshev.dxPhys(b,xMax)
#
#plt.figure(1)
#plt.plot(x,bx1,'b',x,bx2,'m',x,bx3,'r')
#plt.figure(2)
#plt.plot(x,bx2-bx1,'b',x,bx3-bx1,'r')
#plt.show()
#crash

#solution.plotContinuous = False
solution.plot = False
solution.useLongi = not options.noLongi

solution.useSchoofBasal = options.useSchoofBasal

(x1,x2) = extremaBPoly(solution)

print "min, max in b at", x1, x2
if(options.region == 1):
  folder = '%s/retreated_Nx_%i' % (root,Nx)
  xgMins = numpy.linspace(0.7,0.6,5) # Just after b=0
  #xgMaxes = 0.974*numpy.ones((5),float) # Just past the local minimum in b
  xgMaxes = x1*numpy.ones((5),float) # the local minimum
elif(options.region == 2):
  folder = '%s/advanced_Nx_%i' % (root,Nx)
  #xgMins = 1.265*numpy.ones((5),float) # Just before the local maximum in b
  xgMins = x2*numpy.ones((5),float) # The local maximum
  #xgMaxes = numpy.linspace(1.5,1.35,5) 
  xgMaxes = 1.55*numpy.ones((5),float) # Before the cliff gets too steep
else:
  folder = '%s/unstable_Nx_%i' % (root,Nx)
  #xgMins = 0.97*numpy.ones((5),float) # Just before the local maximum in b
  #xgMaxes = 1.27*numpy.ones((5),float) # Before the cliff gets too steep
  xgMins = x1*numpy.ones((5),float) # the local minimum
  xgMaxes = x2*numpy.ones((5),float) # The local maximum


print folder
if not(os.path.exists(folder)):
	os.mkdir(folder)

pIndex = options.pIndex

p = ps[pIndex]
  
for AIndex in range(len(As)):
  A = As[AIndex]
  fileName = '%s/p_%.4f_A_%.4e.da'%(folder,p,A)
  if(os.path.exists(fileName)):
    continue
 
  if(hasattr(solution,'Hprev1')):
    del solution.Hprev1
 
 
  success = solution.findSteadyStateBrentQxg(p, A, xgMins[pIndex], 
            xgMaxes[pIndex], eps_s, tolerance, computeBPoly, computeHGuess)
  if(not success):
    print "Did not find a solution, not writing to a file."
    continue
    
      
  filePointer = open(fileName,'wb')
  solution.Hx.tofile(filePointer)
  solution.H.tofile(filePointer)
  solution.u.tofile(filePointer)
  x = solution.cheb.x*solution.xg
  x.tofile(filePointer)
  solution.xg.tofile(filePointer)
  solution.resStress.tofile(filePointer)
  solution.resBC.tofile(filePointer)
  solution.p.tofile(filePointer)
  solution.A.tofile(filePointer)
  solution.eps_s.tofile(filePointer)
  filePointer.close()

#plt.ioff()
#plt.show()

# vim: ts=2 noet
