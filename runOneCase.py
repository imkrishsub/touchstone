#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 

parser = OptionParser()
parser.add_option("--inFile", type="string", default="none", dest="inFile")
options, args = parser.parse_args()

caseFile = args[0]
caseIndex = int(args[1])
inputFile = options.inFile

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]

lineIndex=3*caseIndex
dir=caseList[lineIndex]
prefix=caseList[lineIndex+1]
commonArgs=caseList[lineIndex+2].split()

print "running case: %s/%s"%(dir,prefix)
if (inputFile != "none"):
 print "reading input from", inputFile

if(not os.path.exists(dir)):
  os.makedirs(dir)
os.chdir(dir)

#subprocess.call(["which", "python"])
#subprocess.call(["whereis", "python"])
#subprocess.call(["ls", "../../code"])

#testArgs="--version"
#print "testing python", testArgs
#subprocess.call(["python", testArgs])
#exit()

#'''
initFile="%s_init"%prefix
looseFile="%s_loose"%prefix
strictFile="%s_strict"%prefix
finalFile="%s_final"%prefix
filePointer="%s.pointer"%prefix

looseTolerance='1e-3'
strictTolerance='1e-6'

commonArgs = ["python", "../../code/mainSheetShelf.py"] \
  + commonArgs \
  + ["--filePointer=%s"%filePointer]

restartFile="none"
if os.path.exists(filePointer):
  lines = [line.rstrip('\n') for line in open(filePointer)]
  restartFile = lines[0]
else:
  if os.path.exists(inputFile):
    restartFile=inputFile
  elif (inputFile != "none"):
    print "exiting: input file", inputFile, "not found."
    exit(1)


if (restartFile == "none"):
  print "init"
  logFile = open("%s.log"%initFile,'w')
  errFile = open("%s.err"%initFile,'w')
  args = commonArgs + ["--outFile=%s.pyda"%initFile, "--eps_s=1e-3", 
    "--maxStep=0", "--maxToleranceInner=1e-3"]
  status = subprocess.call(args, stdout=logFile, stderr=errFile)
  logFile.close()
  errFile.close()
  if status != 0:
    print "init failed! Exiting."
    exit(status)
  restartFile="%s.pyda"%initFile
else:
  print "restarting from file:", restartFile

if (restartFile != "%s.pyda"%strictFile) and ( restartFile != "%s.pyda"%finalFile): 
  print "loose from", restartFile
  logFile = open("%s.log"%looseFile,'w')
  errFile = open("%s.err"%looseFile,'w')
  args = commonArgs + ["--outFile=%s.pyda"%looseFile, "--eps_s=1e-3", 
    "--maxStep=100000", "--maxToleranceInner=1e-3", "--inFile=%s"%restartFile, 
    "--toleranceH=%s"%looseTolerance, "--toleranceXg=%s"%looseTolerance]
  status = subprocess.call(args, stdout=logFile, stderr=errFile)
  logFile.close()
  errFile.close()
  if status != 0:
    print "loose failed! Exiting."
    exit(status)
  restartFile="%s.pyda"%looseFile

if (restartFile == "%s.pyda"%looseFile) or (restartFile == "%s.pyda"%strictFile): 
  print "strict from", restartFile
  logFile = open("%s.log"%strictFile,'w')
  errFile = open("%s.err"%strictFile,'w')
  args = commonArgs + ["--outFile=%s.pyda"%strictFile, "--eps_s=1e-8", 
    "--maxStep=100000", "--maxToleranceInner=1e-5", "--inFile=%s"%restartFile, 
    "--toleranceH=%s"%strictTolerance, "--toleranceXg=%s"%strictTolerance]
  status = subprocess.call(args, stdout=logFile, stderr=errFile)
  logFile.close()
  errFile.close()
  if status != 0:
    print "strict failed! Exiting."
    exit(status)
  restartFile="%s.pyda"%strictFile

if (restartFile == "%s.pyda"%strictFile):
  print "final from", restartFile
  logFile = open("%s.log"%finalFile,'w')
  errFile = open("%s.err"%finalFile,'w')
  args = commonArgs + ["--outFile=%s.pyda"%finalFile, "--eps_s=1e-8", 
    "--maxStep=1000", "--maxToleranceInner=1e-6", "--inFile=%s"%restartFile, 
    "--toleranceH=%s"%strictTolerance, "--toleranceXg=%s"%strictTolerance]
  status = subprocess.call(args, stdout=logFile, stderr=errFile)
  logFile.close()
  errFile.close()
  if status != 0:
    print "final failed! Exiting."
    exit(status)
