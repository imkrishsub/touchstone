#!/usr/bin/python
import numpy
#import math
import matplotlib.pyplot as plt
import os.path
from SteadyStateSolution import SteadyStateSolution
from optparse import OptionParser

def computeBLinear(x,solution):
  shoofx = 750000.
  bx = 778.5*solution.xBar/shoofx/solution.HBar
  b0 = -720./solution.HBar
  eps_b = 1e-8
  absx = numpy.sqrt((bx*x)**2 + eps_b**2)
  b = b0 + absx
  return b

 
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

parser.add_option("--dense", action="store_true", dest="dense", 
                  help="densely sample A?")

parser.add_option("--pIndex", type="int", default=0, dest="pIndex", 
                  help="index into p array (0, 0.25, 0.5, 0.75, 1)")

options, args = parser.parse_args()


Nx = 1025
#eps_s = 1e-3

As = [464.16,
     215.44,
     100.,
     46.416,
     21.544,
     10.,
     4.6416,
     2.1544,
     1.]
As = numpy.array(As)*1e-26 #Pa^(-3) s^(-1)
if(options.dense):
  oldAs = As
  As = numpy.zeros((4*len(As)-3),float)
  As[0::4] = oldAs
  for offset in range(1,4):
    frac = offset/4.0
    As[offset::4] = (1.0 - frac)*oldAs[0:-1]+frac*oldAs[1:] 
  #As = 10.0**numpy.linspace(-24.5,-26.6,41)


eps_s = 1e-8
tolerance = 1e-10

if(options.dense):
  folder = './linear_dense_Ab'
else:
  folder = './linear_Ab_eps_s_1e-8'

print folder
if not(os.path.exists(folder)):
	os.mkdir(folder)

ps = numpy.linspace(0.0,1.0,5)
xgMins = numpy.linspace(1.05,0.85,5) # Before the cliff gets too steep
xgMaxes = 2.0*numpy.ones((5),float) # Just before the local maximum in b

#ps = ps[-1:]
#As = As[-2:]
#xgMax = 0.5*(xgMax+xgMin)

solution = SteadyStateSolution(Nx)
#solution.plotContinuous = False
solution.plot = False

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
            xgMaxes[pIndex], eps_s, tolerance, computeBLinear, computeHGuess)
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
