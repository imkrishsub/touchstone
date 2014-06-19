import numpy

fileName = "allMISMIPCases.txt"

slope_ref = "778.5"
#C_ref = "7.624e6"
lambda_ref = "2.0"
xc = "1.760"
xgInit = "0.8"
eps_s = "1e-10"
m_0 = "0.5"
Ab = "3.1688e-24"
A_ref = "4.6416e-25"

dxs = [ "3.2", "1.6", "0.8", "0.4", "0.2", "0.1", "0.05" ]
Nxs = [ "551", "1101", "2201", "4401", "8801", "17601", "35201" ]
#dts = [ "1e-4", "1e-4", "3e-4", "3e-4", "3e-4", "3e-4", "3e-4" ]
dtInits = [ "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6" ]
ps = [ 0., 0.25, 0.5, 0.75, 1.]
glpStrings = [ "", "--useGLP" ]
glpDirs = [ "nonGLP", "GLP" ]
#As = [ "4.6416e-24", "2.1544e-24", "1.0000e-24", 
#      "4.6416e-25", "2.1544e-25", "1.0000e-25", 
#      "4.6416e-26", "2.1544e-26", "1.0000e-26" ]

Cs = [ "0.316e6", "0.501e6", "0.794e6", "1.258e6", "1.994e6", "3.160e6",
      "5.008e6", "7.624e6", "12.58e6", "19.94e6", "31.60e6", "50.08e6", 
      "79.38e6", "125.8e6" ]

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
tols[2,1,0] = 5e-2 # 0.8 km, p=0.25, nonGLP
tols[2,1,1] = 5e-2 # 0.8 km, p=0.25, GLP
tols[3,1,1] = 5e-2 # 0.4 km, p=0.25, GLP

# p=0.5 at 3.2 and 1.6 km resolution
#goalCFLs[0,2] = 0.01
#goalCFLs[1,2] = 0.01
tols[0,2,1] = 5e-2 # 3.2 km, p=0.5, GLP
tols[1,2,1] = 5e-2 # 1.6 km, p=0.5, GLP

commonArgs = "--folder=. --linearSlope=%s --A=%s --lambda_0=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s"%(slope_ref,A_ref,lambda_ref,eps_s,xgInit,m_0,Ab)
# uncomment the following to include plotting
#commonArgs="%s --plot --plotContinuous"%(commonArgs) 

filePointer = open(fileName, 'w')
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
      exptIndex += 1
      for C in Cs:
        dir = "%s/p_%.2f/res_%s"%(glpDir,p,dx)

        prefix = "C_%s_adv"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s --toleranceInner=%s" \
          %(commonArgs,p,C,dtInit,goalCFL,dx,Nx,xc,glpString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for C in Cs[-2::-1]:
        dir = "%s/p_%.2f/res_%s"%(glpDir,p,dx)

        prefix = "C_%s_ret"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s --toleranceInner=%s" \
          %(commonArgs,p,C,dtInit,goalCFL,dx,Nx,xc,glpString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix


