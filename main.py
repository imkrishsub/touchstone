#!/usr/bin/python
import numpy
from SteadyStateSolver import SteadyStateSolver
from optparse import OptionParser

from pprint import pprint

import os.path

from initH import initH

def computeBLinear(x,solver):
  shoofx = 750000.
  slope = solver.linearSlope/solver.HBar
  pente = 778.5/solver.HBar
  intercept = -720/solver.HBar

  bx = solver.linearSlope*solver.xBar/shoofx/solver.HBar
  b0 = slope*intercept/pente

  eps_b = 1e-8
  absx = numpy.sqrt(x**2 + eps_b**2)
  b = b0 + absx*bx
  
  bx = bx*x/absx
 
  return (b,bx)
  
  
def computeBPoly(x,solver):
  schoofx = 750000.
  xs = x*solver.xBar/schoofx
  b = -(729 - 2184.8*xs**2 + 1031.72*xs**4 - 151.72*xs**6)/solver.HBar  
  bx = -solver.xBar/schoofx*(- 2.0*2184.8*xs + 4.0*1031.72*xs**3 - 6.0*151.72*xs**5)/solver.HBar  
  return (b,bx)

 
#def computeHGuess(x,solver):
#  H0 = 3.0
#  Hf_xg = solver.computeB(x[-1],solver)/(1.0-solver.delta)
#  #print Hf_xg
#  if(Hf_xg <= 0.):
#    raise ValueError("Hf(xg) <= 0. Cannot produce HGuess")
#  H = (Hf_xg - H0)*(x/x[-1])**2.0 + H0
#  
#  return H
  
def readResults(fileName, solver):
  print "reading from: ", fileName
  filePointer = open(fileName,'rb')
  Nx = numpy.fromfile(filePointer, dtype=int, count=1)[0]
  if(Nx != solver.Nx):
    print "Bad Nx in file."
    exit(1)
  
  solver.xg = numpy.fromfile(filePointer, dtype=float, count=1)[0]
  solver.time = numpy.fromfile(filePointer, dtype=float, count=1)[0]
  solver.H = numpy.fromfile(filePointer, dtype=float, count=Nx)
  solver.updateH(solver.H,solver.xg)
  solver.u = numpy.fromfile(filePointer, dtype=float, count=Nx)
  solver.ux = numpy.dot(solver.cheb.Dx/solver.xg,solver.u)
  filePointer.close()

def writeResults(fileName, solver):
  print "writing to: ", fileName
  filePointer = open(fileName,'wb')
  Nx = numpy.array(solver.Nx,int)
  Nx.tofile(filePointer)
  xg = numpy.array(solver.xg)
  xg.tofile(filePointer)
  time = numpy.array(solver.time)
  time.tofile(filePointer)
  solver.H.tofile(filePointer)
  solver.u.tofile(filePointer)
  x = xg*solver.cheb.x
  x.tofile(filePointer)
  solver.Hx.tofile(filePointer)
  solver.ux.tofile(filePointer)
  p = numpy.array(solver.p)
  p.tofile(filePointer)
  A = numpy.array(solver.A)
  A.tofile(filePointer)
  C = numpy.array(solver.C)
  C.tofile(filePointer)
  rho_i = numpy.array(solver.rho_i)
  rho_i.tofile(filePointer)
  a = numpy.array(solver.a)
  a.tofile(filePointer)
  linearSlope = numpy.array(solver.linearSlope)
  linearSlope.tofile(filePointer)
  eps_s = numpy.array(solver.eps_s)
  eps_s.tofile(filePointer)
  filePointer.close()
  filePointer = open("%s/%s"%(options.folder,options.filePointer),'w')
  filePointer.write("%s\n"%options.outFile)
  filePointer.close()


parser = OptionParser()

parser.add_option("--p", type="float", default=0.0, dest="p")
parser.add_option("--A", type="float", default=1e-25, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--rho_i", type="float", default=900.0, dest="rho_i")
parser.add_option("--a", type="float", default=1.0, dest="a")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007
parser.add_option("--lambda_0", type="float", default=2, dest="lambda_0")
parser.add_option("--poly", action="store_true", dest="poly")

parser.add_option("--inFile", type="string", default="none", dest="inFile")
parser.add_option("--folder", type="string", default="results", dest="folder")
parser.add_option("--outFile", type="string", default="results.pyda", dest="outFile")
parser.add_option("--filePointer", type="string", default="default.pointer", dest="filePointer")

parser.add_option("--Nx", type="int", default=1025, dest="Nx")
parser.add_option("--maxSteps", type="int", default=1000, dest="maxSteps")
parser.add_option("--maxInnerSteps", type="int", default=200, dest="maxInnerSteps")
parser.add_option("--maxInitSteps", type="int", default=200, dest="maxInitSteps")
parser.add_option("--stepsPerWrite", type="int", default=10, dest="stepsPerWrite")

parser.add_option("--goalCFL", type="float", default=2.0, dest="goalCFL")
parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")
parser.add_option("--toleranceH", type="float", default=0.1, dest="toleranceH")  # 0.03 m/yr
parser.add_option("--toleranceXg", type="float", default=0.1, dest="toleranceXg")  # 30 m/yr
parser.add_option("--maxToleranceInner", type="float", default=1e-3, dest="maxToleranceInner")
parser.add_option("--initUTolerance", type="float", default=1e-3, dest="initUTolerance")
parser.add_option("--initDt", type="float", default=1e-4, dest="initDt")
parser.add_option("--minXg", type="float", default=0.75, dest="minXg")
parser.add_option("--maxXg", type="float", default=1.2, dest="maxXg")

parser.add_option("--plot", action="store_true", dest="plot")
parser.add_option("--plotContinuous", action="store_true", dest="plotContinuous")

options, args = parser.parse_args()

print "options:"
pprint(options.__dict__)


Nx = options.Nx

eps_s = options.eps_s
toleranceH = options.toleranceH
toleranceXg = options.toleranceXg
maxToleranceInner = options.maxToleranceInner
initUTolerance = options.initUTolerance

stepsPerWrite = options.stepsPerWrite
maxSteps = options.maxSteps
maxInnerSteps = options.maxInnerSteps
maxInitSteps = options.maxInitSteps

goalCFL = options.goalCFL

solver = SteadyStateSolver(Nx)

solver.plotContinuous = options.plotContinuous
solver.plot = options.plot

#solver.useLongi = False

solver.eps_s = eps_s
solver.p = options.p
solver.A = options.A
solver.C = options.C
solver.rho_i = options.rho_i
solver.a = options.a
solver.linearSlope = options.linearSlope
solver.lambda_0 = options.lambda_0

solver.computeNondimConstants()

if(options.poly):
  solver.computeB = computeBPoly
else:
  solver.computeB = computeBLinear

if(not os.path.exists(options.folder)):
  os.mkdir(options.folder)
  
outFile = "%s/%s"%(options.folder,options.outFile)


if(solver.plot):
  import matplotlib.pyplot as plt

if(solver.plot and solver.plotContinuous):
  plt.ion()

if(options.inFile == "none"):
  solver.time = 0
  (HGuess,xgGuess) = initH(solver,options.minXg,options.maxXg)
  
  solver.updateH(HGuess,xgGuess)
  
  solver.u = numpy.zeros(solver.H.shape)
  solver.ux = numpy.zeros(solver.H.shape)
  solver.oldH = solver.H
  solver.oldXg = solver.xg

  innerConverged = False
  print "Computing initial u by iterating on viscosity:"
  for inner in range(maxInitSteps):
    prevU = solver.u
    solver.iterateOnViscosity()
    diffU = numpy.amax(numpy.abs(prevU-solver.u))/numpy.amax(numpy.abs(solver.u))
    if(solver.plot):
      if(solver.plotContinuous):
        plt.draw()
      else:
        plt.show()
    print "iter: ", inner, "diffU:", diffU, "resStress: ", solver.resStress
    if(diffU < initUTolerance):
      innerConverged = True
      break

  if not innerConverged:
    print "Warning: initial velocity did not converge after %i steps!"%maxInitSteps
    
  writeResults(outFile, solver)
  
  # suggest a conservative estimate as the staring time step
  print "Suggested dt:", 0.25*options.goalCFL*(solver.xg/solver.Nx)/numpy.amax(numpy.abs(solver.u))
else:
  if(not os.path.exists(options.inFile)):
    print "File not found:", options.inFile
    exit(1)
  readResults(options.inFile, solver)

if(maxSteps == 0):
  exit()

solver.oldH = solver.H
solver.oldXg = solver.xg

solver.dt = options.initDt

print "initial solver:"
pprint(solver.__dict__)

converged = False

deltaX = numpy.zeros(solver.cheb.x.shape)
deltaX[0:(Nx+1)/2] = solver.cheb.x[1:(Nx+3)/2]-solver.cheb.x[0:(Nx+1)/2]
deltaX[(Nx+1)/2:] = solver.cheb.x[(Nx+1)/2:]-solver.cheb.x[(Nx-1)/2:-1]


toleranceInner = maxToleranceInner
for outer in range(maxSteps):
  solver.oldH = solver.H
  solver.oldXg = solver.xg
  
  innerConverged = False
  
  newH = solver.iterateImplicitTimeStep()
  newXg = solver.newtonStepXg(newH)
  solver.updateH(newH,newXg)

  for inner in range(maxInnerSteps):
    prevU = solver.u
    prevH = solver.H
    solver.iterateOnViscosity()
    newH = solver.iterateImplicitTimeStep()
    newXg = solver.newtonStepXg(newH)

    diffH = numpy.amax(numpy.abs(newH-prevH))/numpy.amax(newH)
    diffXg = numpy.abs(newXg-solver.xg)/newXg
    solver.updateH(newH,newXg)

    diffU = numpy.amax(numpy.abs(prevU-solver.u))/numpy.amax(numpy.abs(solver.u))

    dxg_dt = (newXg- solver.oldXg)/solver.dt
    
    uEff = solver.u - solver.cheb.x*dxg_dt
    cfl = numpy.amax(numpy.abs(uEff/(solver.xg*deltaX)))*solver.dt
    
    if numpy.isnan(cfl):
      print "blew up!"
      exit(1)
      
    print "time: ", solver.time+solver.dt, "iter: ", inner, "cfl:", cfl, "diffs: ", diffH, diffXg, diffU, solver.resStress
    
    if(solver.plot):
      if(solver.plotContinuous):
        plt.draw()
      else:
        plt.show()

    if(diffH < toleranceInner and diffXg < toleranceInner and diffU < toleranceInner):
      innerConverged = True
      break
    
  if not innerConverged:
    print "Error: inner loop did not converge after %i steps!"%maxInnerSteps
    print "Try reducing the goal CFL number."
    exit(1)
    
  dH_dt = (newH-solver.oldH)/solver.dt
  dxg_dt = (newXg-solver.oldXg)/solver.dt
  diffH = numpy.amax(numpy.abs(dH_dt))
  diffXg = numpy.abs(dxg_dt)
  solver.time += solver.dt
  if(diffH < toleranceH and diffXg < toleranceXg):
    converged = True
    break
  # make the inner tolerance stricter as H and xg converge
  # changes in the inner loop should be at most around 10% of changes in the outer loop
  #toleranceInner = min(maxToleranceInner,0.1*solver.dt*max(diffH/numpy.amax(newH),diffXg/newXg))
  toleranceInner = min(maxToleranceInner,solver.dt*max(diffH/numpy.amax(newH),diffXg/newXg))

  scale = max(0.5,min(2.0,goalCFL/cfl))
  solver.dt *= scale

  print "time: ", solver.time, "|dH_dt|_max: ", diffH, "dxg_dt:", dxg_dt, "dt: ", solver.dt, "inner tol.:", toleranceInner
  solver.updateH(newH,newXg)
    
  if(numpy.mod(outer,stepsPerWrite) == stepsPerWrite-1):
    writeResults(outFile, solver)
  
if(converged):
  print "The run converged."
else:
  print "The run did not converge after %i steps."%maxSteps

writeResults(outFile, solver)

if(not converged):
  exit(2)


