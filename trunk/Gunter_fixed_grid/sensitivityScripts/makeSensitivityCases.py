import numpy

fileName = "allSensitivityCases.txt"

slope_ref = "778.5"
C_ref = "7.624e6"
lambda_ref = "2.0"
xc = "2.112"
xgInit = "0.7"
eps_s = "1e-10"
m_0 = "0.5"
A_ref = "2.1544e-25"
Ab = "3.1688e-24"

dxs = [ "3.2", "1.6", "0.8", "0.4", "0.2", "0.1", "0.05" ]
Nxs = [ "661", "1321", "2641", "5281", "10561", "21121", "42241" ]
#dts = [ "1e-4", "1e-4", "3e-4", "3e-4", "3e-4", "3e-4", "3e-4" ]
dtInits = [ "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6", "1e-6" ]
ps = [ 0., 0.25, 0.5, 0.75, 1.]
glpStrings = [ "", "--useGLP" ]
glpDirs = [ "nonGLP", "GLP" ]

#channelStrings = ["", "--useChannelWidth"]
#channelDirs = [ "noButtressing", "withButtressing" ]

#As = [ "4.6416e-24", "2.1544e-24", "1.0000e-24",
#      "4.6416e-25", "2.1544e-25", "1.0000e-25",
#      "4.6416e-26", "2.1544e-26", "1.0000e-26" ]

# We need to start to iterate on the smallest C-value (the least advanced GL).
Cs = [ "3.160e5", "5.008e5", "7.938e5", "1.258e6", "1.994e6",
      "3.160e6", "5.008e6", "7.624e6", "1.258e7", "1.994e7",
      "3.160e7", "7.9384e7", "1.258e8"]
      
 # We need to start to iterate on the steepest slope (the least advanced GL).
slopes = [ "155.70e+2", "778.50e+1", "389.25e+1", "311.40e+1",
          "233.55e+1", "155.70e+1", "778.50e+0", "584.00e+0" ]
      
# We need to start to iterate on the smallest lambda (the least advanced GL).
lambdas = [ "3.125e-2", "6.250e-2", "1.250e-1", "2.500e-1", "5.000e-1",
           "1.000e+0", "2.000e+0", "4.000e+0", "8.000e+0",
           "1.600e+1", "3.200e+1", "6.400e+1", "1.280e+2" ]

# We need to start to iterate on the biggest channel width (the least advanced GL).
Ws = [ "0.5", "1", "2", "3", "4", "5", "10", "20", "100" ]

# We need to start to iterate on the closest calving front (the least advanced GL).
# Remark: when doing this experiment, we need to adjust the values of Nxs!
# xc = [ "1.7", "1.8", "1.9", "2.0", "2.2", "2.4" ]

defaultTol = 1e-3
tols = defaultTol*numpy.ones((len(dxs),len(ps),len(glpStrings)))

# for most cases, we hold the time step roughtly constant
# with resolution (so CFL is proportional to resolution)
# and decrease it by a factor of 2 with increasing p
p0GoalCFLs = numpy.array([ 0.5, 1., 4., 8., 16., 32., 64.])
p1GoalCFLs = p0GoalCFLs
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
#tols[0,1,1] = 5e-2 # 3.2 km, p=0.25, GLP
#tols[1,1,1] = 5e-2 # 1.6 km, p=0.25, GLP
#tols[2,1,0] = 5e-2 # 0.8 km, p=0.25, nonGLP
#tols[2,1,1] = 5e-2 # 0.8 km, p=0.25, GLP
#tols[3,1,1] = 5e-2 # 0.4 km, p=0.25, GLP

# p=0.5 at 3.2 and 1.6 km resolution
#goalCFLs[0,2] = 0.01      
#goalCFLs[1,2] = 0.01
#tols[0,2,1] = 5e-2 # 3.2 km, p=0.5, GLP
#tols[1,2,1] = 5e-2 # 1.6 km, p=0.5, GLP

commonArgs = "--folder=. --A=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s"%(A_ref,eps_s,xgInit,m_0,Ab)
# uncomment the following to include plotting
#commonArgs="%s --plot --plotContinuous"%(commonArgs)

# The first set of loops are without buttressing!
filePointer = open(fileName, 'w')
exptIndex = 0
for resIndex in range(len(dxs)):
  for pIndex in range(len(ps)):
    p = ps[pIndex]
    for glpIndex in range(len(glpStrings)):
      glpDir = glpDirs[glpIndex]
      glpString = glpStrings[glpIndex]
      channelString = ""
#      channelDir = channelDirs[0]
      Nx = Nxs[resIndex]
      dx = dxs[resIndex]
      #dt = dts[resIndex]
      dtInit = dtInits[resIndex]
      goalCFL = "%.2f"%goalCFLs[resIndex,pIndex,glpIndex]
      tol = "%.1e"%tols[resIndex,pIndex,glpIndex]
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
# Iteration over basal stress constant
      for C in Cs:
        dir = "%s/p_%.2f/vary_C/res_%s"%(glpDir,p,dx)      
      
        prefix = "C_%s_adv"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C,slope_ref,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for C in Cs[-2::-1]:
        dir = "%s/p_%.2f/vary_C/res_%s"%(glpDir,p,dx)

        prefix = "C_%s_ret"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C,slope_ref,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

# Iteration over the bed slope topography
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
      for slope in slopes:
        dir = "%s/p_%.2f/vary_Slope/res_%s"%(glpDir,p,dx)

        prefix = "Slope_%s_adv"%slope
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)      
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for slope in slopes[-2::-1]:
        dir = "%s/p_%.2f/vary_Slope/res_%s"%(glpDir,p,dx)

        prefix = "Slope_%s_ret"%slope
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

# Iteration over wavelength
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
      for Lambda in lambdas:
        dir = "%s/p_%.2f/vary_lambda_0/res_%s"%(glpDir,p,dx)

        prefix = "lambda_0_%s_adv"%Lambda
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope_ref,Lambda,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)     
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for Lambda in lambdas[-2::-1]:
        dir = "%s/p_%.2f/vary_lambda_0/res_%s"%(glpDir,p,dx)

        prefix = "lambda_0_%s_ret"%Lambda
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope_ref,Lambda,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

# Iteration over channel width
      channelString = "--useChannelWidth" #channelStrings[1]
      #channelDir = channelDirs[1]
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
      for W in Ws:
        dir = "%s/p_%.2f/vary_W/res_%s"%(glpDir,p,dx)

        prefix = "W_%s_adv"%W
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --W=%s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope_ref,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,W,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix     
      
      for W in Ws[-2::-1]:
        dir = "%s/p_%.2f/vary_W/res_%s"%(glpDir,p,dx)

        prefix = "W_%s_ret"%W
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --W=%s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope_ref,lambda_ref,dtInit,goalCFL,dx,Nx,xc,glpString,channelString,W,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix
     


