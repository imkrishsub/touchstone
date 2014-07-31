import numpy

fileName = "allTransientCases.txt"

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
stepsPerWrite = 2**numpy.arange(0,7)
dtInits = 3e-4/stepsPerWrite
ps = [ 0., 0.25, 0.5, 0.75, 1.]
glpStrings = [ "", "--useGLP" ]
glpDirs = [ "nonGLP", "GLP" ]

# We need to start to iterate on the smallest C-value (the least advanced GL).
Csmall = [ "1.994e6", "7.624e6"]

Clarge = [ "1.994e7", "7.624e6"]
      
defaultTol = 5e-2

commonArgs = "--folder=. --A=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s --maxSteps=1000000 --linearSlope=%s --lambda_0=%s --dtInit=%s --toleranceInner=%s"%(A_ref,eps_s,xgInit,m_0,Ab)
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
      goalCFLvaryC = "%.2f"%goalCFLsvaryC[resIndex,pIndex,glpIndex]
      goalCFLvaryslope = "%.2f"%goalCFLsvaryslope[resIndex,pIndex,glpIndex]
      goalCFLvarylambda = "%.2f"%goalCFLsvarylambda[resIndex,pIndex,glpIndex]
      goalCFLvaryW = "%.2f"%goalCFLsvaryW[resIndex,pIndex,glpIndex]
      tol = "%.1e"%tols[resIndex,pIndex,glpIndex]
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
# Iteration over basal stress constant
      for C in Cs:
        dir = "%s/p_%.2f/vary_C/res_%s"%(glpDir,p,dx)      
      
        prefix = "C_%s_adv"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C,slope_ref,lambda_ref,dtInit,goalCFLvaryC,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for C in Cs[-2::-1]:
        dir = "%s/p_%.2f/vary_C/res_%s"%(glpDir,p,dx)

        prefix = "C_%s_ret"%C
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C,slope_ref,lambda_ref,dtInit,goalCFLvaryC,dx,Nx,xc,glpString,channelString,tol)
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
          %(commonArgs,p,C_ref,slope,lambda_ref,dtInit,goalCFLvaryslope,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)      
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for slope in slopes[-2::-1]:
        dir = "%s/p_%.2f/vary_Slope/res_%s"%(glpDir,p,dx)

        prefix = "Slope_%s_ret"%slope
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope,lambda_ref,dtInit,goalCFLvaryslope,dx,Nx,xc,glpString,channelString,tol)
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
          %(commonArgs,p,C_ref,slope_ref,Lambda,dtInit,goalCFLvarylambda,dx,Nx,xc,glpString,channelString,tol)     
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

      for Lambda in lambdas[-2::-1]:
        dir = "%s/p_%.2f/vary_lambda_0/res_%s"%(glpDir,p,dx)

        prefix = "lambda_0_%s_ret"%Lambda
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --toleranceInner=%s" \
          %(commonArgs,p,C_ref,slope_ref,Lambda,dtInit,goalCFLvarylambda,dx,Nx,xc,glpString,channelString,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix

# Iteration over channel width
      xgInitW = 0.9
      commonArgsW = "--folder=. --A=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s"%(A_ref,eps_s,xgInitW,m_0,Ab)
      channelString = "--useChannel" #channelStrings[1]
      #channelDir = channelDirs[1]
      prevResult = "none"
      print "SENSITIVITY expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
      for W in Ws:
        dir = "%s/p_%.2f/vary_W/res_%s"%(glpDir,p,dx)

        prefix = "W_%s_adv"%W
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --W=%s --toleranceInner=%s" \
          %(commonArgsW,p,C_ref,slope_ref,lambda_ref,dtInit,goalCFLvaryW,dx,Nx,xc,glpString,channelString,W,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix     
      
      for W in Ws[-2::-1]:
        dir = "%s/p_%.2f/vary_W/res_%s"%(glpDir,p,dx)

        prefix = "W_%s_ret"%W
        args = "%s --maxSteps=1000000 --p=%.2f --C=%s --linearSlope=%s --lambda_0=%s --dtInit=%s --goalCFL=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s --W=%s --toleranceInner=%s" \
          %(commonArgsW,p,C_ref,slope_ref,lambda_ref,dtInit,goalCFLvaryW,dx,Nx,xc,glpString,channelString,W,tol)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix
     

