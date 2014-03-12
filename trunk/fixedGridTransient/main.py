#!/usr/bin/python
import numpy
import matplotlib.pyplot as plt
import os.path
from Solver import Solver
from optparse import OptionParser

  
parser = OptionParser()

parser.add_option("--poly", action="store_true", dest="poly", 
                  help="use the polynomial bed")

parser.add_option("--pIndex", type="int", default=0, dest="pIndex", 
                  help="index into p array (0, 0.25, 0.5, 0.75, 1)")

parser.add_option("--AIndex", type="int", default=0, dest="AIndex", 
                  help="index into A array")


parser.add_option("--deltaX", type="float", default=1.6, dest="deltaX", 
                  help="grid spacing in km")

parser.add_option("--deltaT", type="float", default=5, dest="deltaT", 
                  help="time step (years); the code will suggest a deltaT as it goes")

parser.add_option("--adaptStep", action="store_true", dest="adaptStep", 
                  help="adapt the time step over time?")

parser.add_option("--cfl", type="float", default=0.1, dest="cfl", 
                  help="the cfl number at the grounding line, used to adjust the time step")

parser.add_option("--Nx", type="int", default=1250, dest="Nx", 
                  help="Number of H- and u-grid points")

parser.add_option("--stepCount", type="int", default=1000, dest="stepCount", 
                  help="number of time steps to run for")
                  
parser.add_option("--fromLastRun", action="store_true", dest="fromLastRun", 
                  help="restart from  a previous run")

parser.add_option("--initialCondition", type="string", default='', dest="initialCondition", 
                  help="A file containing data to be used as an initial condition (usually the result of a previous run).")

parser.add_option("--continuousPlot", action="store_true", dest="continuousPlot", 
                  help="plot data during simulation")

parser.add_option("--timeCentering", type="float", default=1.0, dest="timeCentering", 
                  help="time centering for thickness evolution -- 1.0 = backward, 0.0 = forward, 0.5 = centered")


options, args = parser.parse_args()


Nx = options.Nx
#eps_s = 1e-3

if(options.poly):
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
else:
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

#if(options.dense):
#  oldAs = As
#  density = 8
#  As = numpy.zeros((density*len(oldAs)-(density-1)),float)
#  As[0::density] = oldAs
#  for offset in range(1,density):
#    frac = offset/float(density)
#    As[offset::density] = (1.0 - frac)*oldAs[0:-1]+frac*oldAs[1:] 


#if(options.dense):
#  denseString = '_dense'
#else:
#  denseString = ''
if(options.poly):
  bedString = 'poly'
else:
  bedString = 'linear'

#folder = './%s%s'%(bedString,denseString)


ps = numpy.linspace(0.0,1.0,5)
pIndex = options.pIndex
p = ps[pIndex]
deltaX = 1e-3*options.deltaX # from km to dimensionless
linearBed = not options.poly

AIndex = options.AIndex

A = As[AIndex]

solver = Solver(Nx,deltaX,p,A,linearBed)
solver.plotContinuous = True
solver.plot = False
solver.transient = True
solver.timeCentering = options.timeCentering

folder = './%s'%(bedString)
if not(os.path.exists(folder)):
	os.mkdir(folder)

folder = '%s/p_%.4f_A_%.4e'%(folder,p,A)
print folder
if not(os.path.exists(folder)):
	os.mkdir(folder)

outputIndex = 1
fileName = '%s/output%i.da'%(folder, outputIndex)
if(options.fromLastRun):
  while(os.path.exists(fileName)):
    outputIndex += 1
    fileName = '%s/output%i.da'%(folder,outputIndex)

if(options.fromLastRun):
  inFileName = '%s/output%i.da'%(folder,outputIndex-1)
else:
  inFileName = options.initialCondition

t0 = 0. 
if(inFileName == ''):
  xg = 1.0
  (H,u) = solver.computeHGuess(xg)
else:
  filePointer = open(inFileName,'rb')
  nxCheck = numpy.fromfile(filePointer,dtype=int,count=1)
  assert(Nx == nxCheck)
  nt = numpy.fromfile(filePointer,dtype=int,count=1)
  H = numpy.fromfile(filePointer,dtype=float,count=Nx)
  u = numpy.fromfile(filePointer,dtype=float,count=Nx)
  if(options.fromLastRun):
    xH = numpy.fromfile(filePointer,dtype=float,count=Nx)
    xu = numpy.fromfile(filePointer,dtype=float,count=Nx)
    oldT = numpy.fromfile(filePointer,dtype=float,count=nt)
    t0 = oldT[-1]
    
  filePointer.close()
  
sPerYr = 365.25*24.*3600. # number of seconds per year
nondimTPerYr = (sPerYr/(solver.xBar/solver.uBar))
deltaT = options.deltaT*nondimTPerYr
#maxT = options.duration*nondimTPerYr
stepCount = options.stepCount #int(numpy.ceil(maxT/deltaT))
maxYears = options.stepCount*options.deltaT
print 'running for %.2f, years'%maxYears
plt.ion()

t = t0 + numpy.array(deltaT*numpy.arange(1,stepCount+1))/nondimTPerYr

xgs = numpy.zeros((stepCount),float)
lambda_g = numpy.zeros((stepCount),float)
deltaH = numpy.zeros((stepCount),float)
newDeltaT = numpy.zeros((stepCount),float)
dtAvg = numpy.zeros((stepCount),float)
glCFL = numpy.zeros((stepCount),float)
for tIndex in range(stepCount):
  print 'step:', tIndex, '/', (stepCount-1)
  HPrev = H
  (H,u) = solver.step(H,u,H,deltaT)
  if options.adaptStep and tIndex > 1:
      t[tIndex] = t[tIndex-1] + deltaT/nondimTPerYr
      t[-1] = t[tIndex]+deltaT/nondimTPerYr*(stepCount-1-tIndex)
  print 'time:', t[tIndex]
  glCFL[tIndex] = numpy.amax(numpy.abs(u[solver.glIndices]))*deltaT/solver.deltaX
  print 'cfl at xg:', glCFL[tIndex]
  newDeltaT[tIndex] = deltaT*options.cfl/glCFL[tIndex]
  
  lambda_g[tIndex] = solver.lambda_g[0]
  print 'lambda_g: ',lambda_g[tIndex]
  xgPrev = 1000.*(solver.xH[solver.glIndicesPrev[0]]+solver.deltaX*solver.lambda_gPrev[0])
  xgs[tIndex] = 1000.*(solver.xH[solver.glIndices[0]]+solver.deltaX*solver.lambda_g[0])
  print 'xg (km): ',xgs[tIndex]
  deltaXg = xgs[tIndex]-xgPrev
  print 'deltaXg (km): ', deltaXg
  #deltaT_xg = options.cfl*deltaT*deltaX/numpy.abs(deltaXg/1000.)
  #newDeltaT[tIndex] = numpy.minimum(deltaT_xg,newDeltaT[tIndex])
    
  window = numpy.minimum(20,tIndex+1)
  dtList = newDeltaT[tIndex-window+1:tIndex+1]
  windowMin = numpy.amin(dtList)
  weights = (windowMin/dtList)**3 #weighted to favor the minimum
  dtAvg[tIndex] = numpy.sum(dtList*weights)/numpy.sum(weights)
  print 'suggested deltaT (years):', dtAvg[tIndex]/nondimTPerYr
  
  if options.adaptStep:
    deltaT = dtAvg[tIndex]

  deltaH[tIndex] = numpy.amax(numpy.abs(H-HPrev))*1000.
  print 'max deltaH (m):', deltaH[tIndex]
    
  if tIndex == 0:
    H0 = H
    u0 = u
    (s0,sx) = solver.computeSx(H,solver.floatingMaskH,
        solver.floatingMaskU,solver.glIndices,solver.lambda_g)
    floatingMaskH0 = solver.floatingMaskH
    
  

  if(options.continuousPlot or tIndex == stepCount-1):
    fig = plt.figure(1)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    plt.plot(solver.xu,u,'b', solver.xu,u0, 'r')        
    plt.axis('tight')
    plt.title('p=%.2f'%p)

    fig = plt.figure(2)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    plt.plot(solver.xu,solver.longi, 'b', solver.xu, solver.basal, 'r', solver.xu, solver.driving, 'g')
    plt.axis('tight')
    plt.title('p=%.2f'%p)
    
    (s,sx) = solver.computeSx(H,solver.floatingMaskH,
        solver.floatingMaskU,solver.glIndices,solver.lambda_g)
    fig = plt.figure(3)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    plt.plot(solver.xH,s,'b', solver.xH, s0, 'r', solver.xH, -solver.b, 'k',
             solver.xH[solver.floatingMaskH],s[solver.floatingMaskH]-H[solver.floatingMaskH], 'b',
             solver.xH[floatingMaskH0],s0[floatingMaskH0]-H0[floatingMaskH0], 'r')
    plt.axis('tight')
    plt.title('p=%.2f'%p)
    
    fig=plt.figure(4)
    for index in range(len(fig.axes)):
      fig.axes[index].cla()
    ax1 = fig.add_subplot(111)
    xgPlot = xgs[:tIndex+1]
    tPlot = t[:tIndex+1]
    ax1.plot(tPlot, deltaH[:tIndex+1], 'r')
    ax1.set_xlim([t[0],t[-1]])
    ax1.set_xlabel('t (years)')
    ax1.set_ylabel('deltaH (m)',color='r')
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
        
    ax2 = ax1.twinx()
    ax2.plot(tPlot, xgPlot,'b')
    ax2.set_ylabel('xg (km)',color='b')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
    plt.title('p=%.2f'%p)
  
    fig = plt.figure(5)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    #plt.plot(tPlot,newDeltaT[:tIndex+1],'b',tPlot,dtAvg[:tIndex+1],'k')
    plt.plot(tPlot,glCFL[:tIndex+1],'b')
    plt.xlim([t[0],t[-1]])
    plt.ylabel('cfl at gl')
    plt.title('p=%.2f'%p)

    fig = plt.figure(6)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    #plt.plot(tPlot,newDeltaT[:tIndex+1],'b',tPlot,dtAvg[:tIndex+1],'k')
    plt.plot(tPlot,lambda_g[:tIndex+1],'b')
    plt.xlim([t[0],t[-1]])
    plt.ylabel('lambda_g')
    plt.title('p=%.2f'%p)
  
    xg = xgs[tIndex]
    x = 1000.*solver.xu-xg
    mask = numpy.logical_and(x <= 10, x >= -40)
    x = x[mask]
    fig = plt.figure(7)
    if(len(fig.axes) > 0):
      fig.axes[0].cla()
    plt.plot(x,solver.longi[mask], 'b', x, solver.basal[mask], 'r', x, solver.driving[mask], 'g')
    plt.axis('tight')
    plt.title('p=%.2f'%p)

    plt.draw()
  
        
filePointer = open(fileName,'wb')
numpy.array(Nx).tofile(filePointer)
numpy.array(len(t)).tofile(filePointer)
H.tofile(filePointer)
u.tofile(filePointer)
solver.xH.tofile(filePointer)
solver.xu.tofile(filePointer)
t.tofile(filePointer)
xgs.tofile(filePointer)
solver.p.tofile(filePointer)
solver.A.tofile(filePointer)
filePointer.close()
  
plt.ioff()

plt.show()


