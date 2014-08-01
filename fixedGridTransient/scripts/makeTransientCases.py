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
Cschoof = "7.624e6"

meltRates = ["1","10","50","100","150"]
      
defaultTol = 5e-2

commonArgs = "--folder=. --writeToSeparateFile --initFromCheb --fixedTimeStep --A=%s --eps_s=%s --xgInit=%s --m_0=%s --Ab=%s --maxSteps=1000000 --linearSlope=%s --lambda_0=%s --toleranceInner=%s"%(A_ref,eps_s,xgInit,m_0,Ab,slope_ref,lambda_ref,defaultTol)


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
#      goalCFLvaryC = "%.2f"%goalCFLsvaryC[resIndex,pIndex,glpIndex]
#      tol = "%.1e"%tols[resIndex,pIndex,glpIndex]

#      prevResult = "none"
      prevResult = "C_%s_cheby.pyda"%Cschoof
      print "TRANSIENT expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1
     
      dir = "%s/p_%.2f/transient_C/res_%s"%(glpDir,p,dx)      
      C = Clarge[0]
      prefix = "C_%s_adv"%C
      args = "%s --p=%.2f --C=%s --dtInit=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s"%(commonArgs,p,C,dtInit,dx,Nx,xc,glpString,channelString)
      filePointer.write("%s\n"%dir)
      filePointer.write("%s\n"%prefix)
      filePointer.write("%s\n"%prevResult)
      filePointer.write("%s\n"%args)
      prevResult="%s_final.pyda"%prefix


      prevResult = "C_%s_cheby.pyda"%Csmall[0]
      print "TRANSIENT expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1

      dir = "%s/p_%.2f/transient_C/res_%s"%(glpDir,p,dx)
      C = Cschoof
      prefix = "C_%s_adv"%C
      args = "%s --p=%.2f --C=%s --dtInit=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s"%(commonArgs,p,C,dtInit,dx,Nx,xc,glpString,channelString)
      filePointer.write("%s\n"%dir)
      filePointer.write("%s\n"%prefix)
      filePointer.write("%s\n"%prevResult)
      filePointer.write("%s\n"%args)
      prevResult="%s_final.pyda"%prefix


      prevResult = "C_%s_cheby.pyda"%Clarge[0]
      print "TRANSIENT expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1

      dir = "%s/p_%.2f/transient_C/res_%s"%(glpDir,p,dx)
      C = Cschoof
      prefix = "C_%s_ret"%C
      args = "%s --p=%.2f --C=%s --dtInit=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s"%(commonArgs,p,C,dtInit,dx,Nx,xc,glpString,channelString)
      filePointer.write("%s\n"%dir)
      filePointer.write("%s\n"%prefix)
      filePointer.write("%s\n"%prevResult)
      filePointer.write("%s\n"%args)
      prevResult="%s_final.pyda"%prefix


      prevResult = "C_%s_cheby.pyda"%Cschoof
      print "TRANSIENT expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
      exptIndex += 1

      dir = "%s/p_%.2f/transient_C/res_%s"%(glpDir,p,dx)
      C = Csmall[0]
      prefix = "C_%s_ret"%C
      args = "%s --p=%.2f --C=%s --dtInit=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s"%(commonArgs,p,C,dtInit,dx,Nx,xc,glpString,channelString) 
      filePointer.write("%s\n"%dir)
      filePointer.write("%s\n"%prefix)
      filePointer.write("%s\n"%prevResult)
      filePointer.write("%s\n"%args)
      prevResult="%s_final.pyda"%prefix


###################################################
      for meltRate in meltRates:
# Iteration over meltrates
        prevResult = "C_%s_cheby.pyda"%Cschoof
        print "TRANSIENT expt %i: %6s, p=%.2f, res=%s"%(exptIndex,glpDir,p,dx)
        exptIndex += 1

        dir = "%s/p_%.2f/transient_meltRate/res_%s"%(glpDir,p,dx)

        prefix = "meltRate_%s"%meltRate
        args = "%s --p=%.2f --C=%s --dtInit=%s --meltRate=%s --deltaX=%se-3 --Nx=%s --xc=%s %s %s"%(commonArgs,p,C_ref,dtInit,meltRate,dx,Nx,xc,glpString,channelString)
        filePointer.write("%s\n"%dir)
        filePointer.write("%s\n"%prefix)
        filePointer.write("%s\n"%prevResult)
        filePointer.write("%s\n"%args)
        prevResult="%s_final.pyda"%prefix



