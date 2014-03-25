#!/usr/bin/python
import numpy
from Solver2 import Solver
from optparse import OptionParser

from pprint import pprint

import os.path
import numpy.linalg


  
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
#  solver.ux = numpy.dot(solver.Dxu,solver.u)
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
  xH = solver.xH
  xH.tofile(filePointer)
  xu = solver.xu
  xu.tofile(filePointer)
#  solver.Hx.tofile(filePointer)
#  solver.ux.tofile(filePointer)
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

parser.add_option("--p", type="float", default=0.25, dest="p")
parser.add_option("--A", type="float", default=1e-24, dest="A")
parser.add_option("--C", type="float", default=7.624e6, dest="C")
parser.add_option("--W", type="float", default=10, dest="W")
parser.add_option("--rho_i", type="float", default=900.0, dest="rho_i")
parser.add_option("--a", type="float", default=1.0, dest="a")
parser.add_option("--linearSlope", type="float", default=778.5, dest="linearSlope") #drop in m per 750 km, as in Schoof 2007
parser.add_option("--lambda_0", type="float", default=2, dest="lambda_0")
parser.add_option("--poly", action="store_true", dest="poly")
parser.add_option("--useGLP", action="store_false", dest="useGLP")

parser.add_option("--inFile", type="string", default="none", dest="inFile")
parser.add_option("--folder", type="string", default="results", dest="folder")
parser.add_option("--outFile", type="string", default="results.pyda", dest="outFile")
parser.add_option("--filePointer", type="string", default="default.pointer", dest="filePointer")

parser.add_option("--Nx", type="int", default=1321, dest="Nx")
parser.add_option("--xc", type="float", default=2.112, dest="xc")
parser.add_option("--deltaX", type="float", default=1.6, dest="deltaX")
parser.add_option("--maxSteps", type="int", default=2000, dest="maxSteps")
parser.add_option("--maxInnerSteps", type="int", default=200, dest="maxInnerSteps")
parser.add_option("--maxSteps", type="int", default=200, dest="maxInitSteps")
parser.add_option("--stepsPerWrite", type="int", default=10, dest="stepsPerWrite")

parser.add_option("--goalCFL", type="float", default=2.0, dest="goalCFL")
parser.add_option("--eps_s", type="float", default=1e-8, dest="eps_s")
parser.add_option("--toleranceH", type="float", default=0.1, dest="toleranceH")  # 0.03 m/yr
parser.add_option("--toleranceXg", type="float", default=0.1, dest="toleranceXg")  # 30 m/yr
parser.add_option("--maxToleranceInner", type="float", default=3e-3, dest="maxToleranceInner")
parser.add_option("--initUTolerance", type="float", default=1e-3, dest="initUTolerance")
parser.add_option("--dt", type="float", default=3e-4, dest="dt")
parser.add_option("--XgInit", type="float", default=1.62, dest="XgInit")

parser.add_option("--plot", action="store_true", dest="plot")
parser.add_option("--plotContinuous", action="store_true", dest="plotContinuous")

options, args = parser.parse_args()

print "options:"
pprint(options.__dict__)

xc = options.xc # calving front on u-grid
deltaX = 1e-3*options.deltaX
Nx = options.Nx
A = options.A
p = options.p
XgInit = options.XgInit

eps_s = options.eps_s
toleranceH = options.toleranceH
toleranceXg = options.toleranceXg
maxToleranceInner = options.maxToleranceInner
initUTolerance = options.initUTolerance

stepsPerWrite = options.stepsPerWrite
maxSteps = options.maxSteps
maxInitSteps = options.maxInitSteps

goalCFL = options.goalCFL

linearBed = not options.poly
linearSlope = options.linearSlope
useGLP = options.useGLP

solver = Solver(Nx,deltaX,xc,p,A,linearBed,linearSlope,useGLP)

solver.plotContinuous = options.plotContinuous
solver.plot = options.plot
solver.maxPicardIter = options.maxInnersteps
#solver.useLongi = False

solver.eps_s = eps_s
solver.p = options.p
solver.A = options.A
solver.C = options.C
solver.rho_i = options.rho_i
solver.a = options.a
solver.linearSlope = options.linearSlope
solver.lambda_0 = options.lambda_0
#solver.xg = options.XgInit


if(not os.path.exists(options.folder)):
  os.mkdir(options.folder)
  
outFile = "%s/%s"%(options.folder,options.outFile)


if(solver.plot):
  import matplotlib.pyplot as plt

if(solver.plot and solver.plotContinuous):
  plt.ion()

if(options.inFile == "none" or options.inFile == "zero"):
  solver.time = 0
  if(options.inFile == "none"):
    (HGuess,uGuess) = solver.computeHGuess(XgInit)
  
  
  solver.H = HGuess
  solver.u = uGuess  
  
# solver.u = numpy.zeros(solver.H.shape)
# solver.ux = numpy.zeros(solver.H.shape)
  solver.oldH = solver.H
  solver.oldXg = solver.xg

 # Check (reinitializing) uGuess using Picard iteration 
  innerConverged = False
  print "Computing initial u by iterating on viscosity:"
  for inner in range(maxInitSteps):
    prevU = solver.u
    solver.iteratePicardTau(solver.H, solver.u)
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
  print "Suggested dt:", 0.25*options.goalCFL*(deltaX)/numpy.amax(numpy.abs(solver.u))
else:
  if(not os.path.exists(options.inFile)):
    print "File not found:", options.inFile
    exit(1)
  readResults(options.inFile, solver)

if(maxSteps == 0):
  exit()

solver.oldH = solver.H    # meaning previous time step
solver.oldXg = solver.xg  # meaning previous time step

dt = options.dt

print "initial solver:"
pprint(solver.__dict__)

converged = False

toleranceInner = maxToleranceInner

#tol_out_H = 1e-1*solver.dt*solver.tBar*solver.aBar/solver.HBar
#tol_out_xg = 1e-3*solver.dt*solver.tBar*solver.aBar/solver.xBar
tol_out_dH_dt = 1e-1/solver.sPerY/solver.aBar # tolerance in non-dimensional unit
tol_out_dxg_dt = 1e-3/solver.sPerY/solver.uBar



for outer in range(maxSteps):
#  solver.oldH = solver.H
#  solver.oldXg = solver.xg
#  solver.oldU = solver.u
  
  innerConverged = False
  
#  newH = solver.solveCont(solver.H, solver.u, HPrev, solver.dt, firstIter)
  solver.xg = solver.updateXg(solver.H,solver.Hf)

  prevU = solver.u
  prevH = solver.H
  prevXg = solver.xg
  
  (Hkp1,ukp1) = solver.step(solver.H, solver.u, prevH, dt)
  
  solver.xg = solver.updateXg(Hkp1,solver.Hf)
  
  #    diffH = numpy.amax(numpy.abs(newH-prevH))/numpy.amax(newH)
  #    diffXg = numpy.abs(newXg-solver.xg)/newXg
  #    diffU = numpy.amax(numpy.abs(prevU-solver.u))/numpy.amax(numpy.abs(solver.u))
  
  dH_dt = numpy.linalg.norm(Hkp1-prevH)/dt
  dxg_dt = numpy.abs(solver.xg-prevXg)/dt
  solver.u = ukp1
  solver.H = Hkp1
  
  #    dxg_dt = (newXg- solver.oldXg)/solver.dt
  cfl = numpy.amax(numpy.abs(solver.u/deltaX))*dt
  print "dH_dt = ",dH_dt, "dxg_dt = ",dxg_dt, "CFL = ",cfl
  
  if numpy.isnan(cfl):
    print "blew up!"
    exit(1)
      
#    print "time: ", solver.time+solver.dt, "iter: ", inner, "cfl:", cfl, "diffs: ", diffH, diffXg, diffU, solver.resStress
        
  if not solver.innerConverged:
    print "Error: inner loop did not converge after %i steps!"%solver.maxPicardIter
    print "Try reducing the goal CFL number."
#    print "diffU=,",diffU
    exit(1)
  
#  newH = solver.solveCont(solver.H, solver.u, HPrev, solver.dt, firstIter)
#  solver.H = newH  
#  newXg = solver.updateXg(newH,solver.Hf)
#  solver.xg = newXg  
  
#  dH_dt = (newH-solver.oldH)/solver.dt
#  prev_dH_dt = dH_dt
#  prev_du_dt = (solver.u-solver.oldU)/solver.dt
#  dxg_dt = (newXg-solver.oldXg)/solver.dt
  #prev_dxg_dt = dxg_dt
  #print "prev_dxg_dt:", prev_dxg_dt
#  diffH = numpy.amax(numpy.abs(dH_dt))
#  diffXg = numpy.abs(dxg_dt)
  solver.time += dt
  if(dH_dt < tol_out_dH_dt and dxg_dt < tol_out_dxg_dt):
    converged = True
    break
  # make the inner tolerance stricter as H and xg converge
  # changes in the inner loop should be at most the size of changes in the outer loop,
  # once we have reached the goal CFL
#  maxChange = max(diffH/numpy.amax(newH),diffXg/newXg)
#  toleranceInner = min(maxToleranceInner,solver.dt*goalCFL/cfl*maxChange)

#  scale = max(0.5,min(2.0,goalCFL/cfl))
#  solver.dt *= scale

#  print "time: ", solver.time, "|dH_dt|_max: ", diffH, "dxg_dt:", dxg_dt, "dt: ", solver.dt, "inner tol.:", toleranceInner
#  solver.updateH(newH,newXg)
#  solver.H = solver.solveCont(newH, solver.u, solver.H, solver.dt, firstIter)
#  solver.xg = newXg
  
  
  if(numpy.mod(outer,stepsPerWrite) == stepsPerWrite-1):
    writeResults(outFile, solver)
  
if(converged):
  print "The run converged."
  print "xg =",solver.xg
else:
  print "The run did not converge after %i steps."%maxSteps
  print "xg =",solver.xg
  
writeResults(outFile, solver)

if(not converged):
  exit(2)


