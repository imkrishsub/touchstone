#!/usr/bin/python
import numpy
#import math
import matplotlib.pyplot as plt
import os.path
from SteadyStateSolution import SteadyStateSolution
from optparse import OptionParser

def computeBLinear(x,solution):
      shoofx = 750000.
      slope = solution.linearSlope/solution.HBar
      pente = 778.5/solution.HBar
      intercept = -720/solution.HBar

      bx = solution.linearSlope*solution.xBar/shoofx/solution.HBar
      b0 = slope*intercept/pente
      eps_b = 1e-8
      absx = numpy.sqrt(x**2 + eps_b**2)
      b = b0 + absx*bx

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

parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=2.1544e-25, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--lambda_0", type="float", default=2.0, dest="lambda_0")
parser.add_option("--slope", type="float", default=778.5, dest="slope")
parser.add_option("--outFile", type="string", default="results/output.pyda", dest="outFile")

options, args = parser.parse_args()


Nx = 1025
eps_s = 1e-8
tolerance = 1e-10

folder = os.path.dirname(options.outFile)

print folder
if not(os.path.exists(folder)):
	os.mkdir(folder)

ps = numpy.linspace(0.0,1.0,5)
xgMin = 0.7
xgMax = 3.0

solution = SteadyStateSolution(Nx)
A = options.A
p = options.p
solution.C = options.C
solution.lambda_0 = options.lambda_0
solution.linearSlope = options.slope
#solution.plotContinuous = False
solution.plot = False

pIndex = options.pIndex

fileName = options.outFile 
success = solution.findSteadyStateBrentQxg(p, A, xgMin,
           xgMax, eps_s, tolerance, computeBLinear, computeHGuess)
if(not success):
  print "Did not find a solution, not writing to a file."
  exit()
    
      
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

