import numpy

fileName = "allMISMIPCases.txt"

slope_ref = "778.5"
C_ref = "7.624e6"
lambda_ref = "2.0"
xc = "1.760"


dxs = [ "3.2", "1.6", "0.8", "0.4", "0.2", "0.1", "0.05" ]
Nxs = [ "551", "1101", "2201", "4401", "8801", "17601", "35201" ]
#dts = [ "1e-4", "1e-4", "3e-4", "3e-4", "3e-4", "3e-4", "3e-4" ]
dtInits = [ "1e-5", "1e-5", "1e-5", "1e-5", "1e-5", "1e-5", "1e-5" ]
p0GoalCFLs = numpy.array([ 0.25, 0.5, 4., 8., 16., 32., 64.])
p1GoalCFLs = 0.5*p0GoalCFLs

As = [ "4.6416e-24", "2.1544e-24", "1.0000e-24", 
      "4.6416e-25", "2.1544e-25", "1.0000e-25", 
      "4.6416e-26", "2.1544e-26", "1.0000e-26" ]
#ps = [ "0.00", "0.25", "0.50", "0.75", "1.00" ]
ps = [ 0., 0.25, 0.5, 0.75, 1.]
glpStrings = [ "", "--useGLP" ]
glpDirs = [ "nonGLP", "GLP" ]
commonArgs = "--folder=. --linearSlope=%s --C=%s --lambda_0=%s --eps_s=1e-8"%(slope_ref,C_ref,lambda_ref)
# uncomment the following to include plotting
#commonArgs="%s --plot --plotContinuous"%(commonArgs) 

filePointer = open(fileName, 'w')
for resIndex in range(len(dxs)):
  for p in ps:
    for glpIndex in range(len(glpStrings)):
      glpDir = glpDirs[glpIndex]
      glpString = glpStrings[glpIndex]
      Nx = Nxs[resIndex]
      dx = dxs[resIndex]
      #dt = dts[resIndex]
      dtInit = dtInits[resIndex]
      goalCFL = "%.2f"%((1.-p)*p0GoalCFLs[resIndex]+p*p1GoalCFLs[resIndex])
      prevResult = "none"
      for A in As:
        dir = "%s/p_%.2f/res_%s"%(glpDir,p,dx)

        prefix = "A_%s_adv"%A
        args = "%s --maxSteps=1000000 --p=%.2f --A=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s" \
          %(commonArgs,p,A,dtInit,goalCFL,dx,Nx,xc,glpString)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for A in As[-2::-1]:
        dir = "%s/p_%s/res_%s"%(glpDir,p,dx)

        prefix = "A_%s_ret"%A
        args = "%s --maxSteps=1000000 --p=%.2f --A=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s" \
          %(commonArgs,p,A,dtInit,goalCFL,dx,Nx,xc,glpString)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix


