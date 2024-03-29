import numpy

fileName = "allPolyMISMIPCases.txt"

C_ref = "7.624e6"
lambda_ref = "2.0"
xc = "1.504"
xgInit = "0.6"
eps_s = "1e-10"
m_0 = "0.5"
Ab = "3.1688e-24"

dxs = [ "3.2", "1.6", "0.8", "0.4", "0.2", "0.1", "0.05" ]
Nxs = [ "471", "941", "1881", "3761", "7521", "15041", "30081" ]
#dts = [ "1e-4", "1e-4", "3e-4", "3e-4", "3e-4", "3e-4", "3e-4" ]
dtInits = [ "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6" ]
ps = [ 0., 0.25, 0.5, 0.75, 1.]
glpStrings = [ "", "--useGLP" ]
glpDirs = [ "nonGLP", "GLP" ]

######################### Gunter put A values here #####################
A_p0 = [ "3.0000e-25", "2.8231e-25", "2.6566e-25",
       "2.5000e-25", "2.3210e-25", "2.1544e-25", 
       "2.0000e-25", "1.8171e-25", "1.6510e-25",
       "1.5000e-25", "1.3103e-25", "1.1447e-25",
       "1.0000e-25", "7.9370e-26", "6.2996e-26",
       "5.0000e-26", "3.9685e-26", "3.1498e-26",
       "2.5000e-26" ]
A_p1 = [ "3.0000e-25", "2.8231e-25", "2.6566e-25", 
      "2.5000e-25", "2.3210e-25", "2.1544e-25", 
      "2.0000e-25", "1.8171e-25", "1.6510e-25",
      "1.5000e-25", "1.3103e-25", "1.1447e-25",
      "1.0000e-25", "7.9370e-26", "6.2996e-26",
      "5.0000e-26", "3.9685e-26", "3.1498e-26",
      "2.5000e-26", "2.3208e-26", "2.1544e-26",
      "2.0000e-26", "1.8171e-26", "1.6510e-26",
      "1.5000e-26", "1.3104e-26", "1.1447e-26",
      "1.0000e-26", "7.9279e-27", "6.2996e-27",
      "5.0000e-27", "3.9685e-27", "3.1498e-27",
      "2.5000e-27"]

As = []
As.append(A_p0)
As.append(A_p0)
As.append(A_p0)
As.append(A_p1)
As.append(A_p1)

defaultTol = 1e-3
tols = defaultTol*numpy.ones((len(dxs),len(ps),len(glpStrings)))

# for most cases, we hold the time step roughtly constant
# with resolution (so CFL is proportional to resolution) 
# and decrease it by a factor of 2 with increasing p
p0GoalCFLs = numpy.array([ 0.5, 1., 4., 8., 16., 32., 64.])
p1GoalCFLs = 0.5*p0GoalCFLs
goalCFLs = numpy.zeros((len(dxs),len(ps),len(glpStrings)))
for resIndex in range(len(dxs)):
  for pIndex in range(len(ps)):
    for glpIndex in range(len(glpStrings)):
      p = ps[pIndex]
      goalCFLs[resIndex,pIndex] = \
        ((1.-p)*p0GoalCFLs[resIndex]+p*p1GoalCFLs[resIndex])

# a few cases are very stubborn and require special treatment
# p=0.25 at 3.2 and 1.6 km resolution
#goalCFLs[0,1] = 0.0001
#goalCFLs[1,1] = 0.0001
tols[0,1,1] = 5e-2 # 3.2 km, p=0.25, GLP
tols[1,1,1] = 5e-2 # 1.6 km, p=0.25, GLP

# p=0.5 at 3.2 and 1.6 km resolution
#goalCFLs[0,2] = 0.01
#goalCFLs[1,2] = 0.01
tols[0,2,1] = 5e-2 # 3.2 km, p=0.5, GLP
tols[1,2,1] = 5e-2 # 1.6 km, p=0.5, GLP

commonArgs = "--poly --folder=. --C=%s --lambda_0=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s"%(C_ref,lambda_ref,eps_s,xgInit,m_0,Ab)
# uncomment the following to include plotting
#commonArgs="%s --plot --plotContinuous"%(commonArgs) 

caseCount = len(dxs)*len(ps)*len(glpStrings)
filePointer = open(fileName, 'w')
filePointer.write("%i\n"%caseCount)

exptIndex = 0
for resIndex in range(len(dxs)):
  for pIndex in range(len(ps)):
    p = ps[pIndex]
    for glpIndex in range(len(glpStrings)):
      glpDir = glpDirs[glpIndex]
      glpString = glpStrings[glpIndex]
      Nx = Nxs[resIndex]
      dx = dxs[resIndex]
      #dt = dts[resIndex]
      dtInit = dtInits[resIndex]
      goalCFL = "%.2f"%goalCFLs[resIndex,pIndex,glpIndex]
      tol = "%.1e"%tols[resIndex,pIndex,glpIndex]
      prevResult = "none"
      print "MISMIP expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      filePointer.write("%i\n"%(2*len(As[pIndex])-1))
      exptIndex += 1
      for A in As[pIndex]:
        dir = "%s/p_%.2f/res_%s"%(glpDir,p,dx)

        prefix = "A_%s_adv"%A
        args = "%s --maxSteps=1000000 --p=%.2f --A=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s --toleranceInner=%s" \
          %(commonArgs,p,A,dtInit,goalCFL,dx,Nx,xc,glpString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for A in As[pIndex][-2::-1]:
        dir = "%s/p_%.2f/res_%s"%(glpDir,p,dx)

        prefix = "A_%s_ret"%A
        args = "%s --maxSteps=1000000 --p=%.2f --A=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s --toleranceInner=%s" \
          %(commonArgs,p,A,dtInit,goalCFL,dx,Nx,xc,glpString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix


