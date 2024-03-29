#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 

parser = OptionParser()
#parser.add_option("--inFile", type="string", default="none", dest="inFile")
options, args = parser.parse_args()

caseFile = args[0]
caseIndex = int(args[1])
lineIndex = int(args[2])
#inputFile = options.inFile

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]

dir=caseList[lineIndex]
prefix=caseList[lineIndex+1]
inputFile=caseList[lineIndex+2]
commonArgs=caseList[lineIndex+3].split()

print "running case: %s/%s"%(dir,prefix)
if (inputFile != "none"):
 print "reading input from", inputFile

if(not os.path.exists(dir)):
  os.makedirs(dir)

oldPath = os.getcwd()

os.chdir(dir)

inProgressFile="%s_inProgress"%prefix
finalFile="%s_final"%prefix
filePointer="%s.pointer"%prefix

commonArgs = ["python", "%s/code/main.py"%oldPath] \
  + commonArgs \
  + ["--folder=.", "--filePointer=%s"%filePointer]

restartFile="none"
if os.path.exists(filePointer):
  lines = [line.rstrip('\n') for line in open(filePointer)]
  if len(lines) > 0 and os.path.exists(lines[0]) and (os.stat(lines[0]).st_size > 0):
    restartFile=lines[0]
if(restartFile == "none"):
  if os.path.exists(inputFile):
    restartFile=inputFile
  elif (inputFile != "none"):
    print "exiting: input file", inputFile, "not found."
    exit(1)

if (restartFile == "%s.pyda"%finalFile):
  print "already done"
  exit(0)

if (restartFile == "none"):
  print "init"
else:
  print "restarting from file:", restartFile

logFile = open("%s.log"%inProgressFile,'w')
errFile = open("%s.err"%inProgressFile,'w')
args = commonArgs + ["--outFile=%s.pyda"%inProgressFile, "--inFile=%s"%restartFile]
status = subprocess.call(args, stdout=logFile, stderr=errFile)
logFile.close()
errFile.close()
if status == 40:
  print "failed to read input file. Trying again with %s."%inputFile

  logFile = open("%s.log"%inProgressFile,'w')
  errFile = open("%s.err"%inProgressFile,'w')
  args = commonArgs + ["--outFile=%s.pyda"%inProgressFile, "--inFile=%s"%inputFile]
  status = subprocess.call(args, stdout=logFile, stderr=errFile)
  logFile.close()
  errFile.close()

if status != 0:
  print "run failed! Exiting."
  exit(status)
#if status == 2:
#  print "warning: run did not converge.  Consider increasing maxSteps."
restartFile="%s.pyda"%inProgressFile

print "final"
logFile = open("%s.log"%finalFile,'w')
errFile = open("%s.err"%finalFile,'w')
args = commonArgs + ["--outFile=%s.pyda"%finalFile, "--inFile=%s"%restartFile]
status = subprocess.call(args, stdout=logFile, stderr=errFile)
logFile.close()
errFile.close()
if status != 0:
  print "final failed! Exiting."
  exit(status)
#if status == 2:
#  print "warning: run did not converge.  Consider increasing maxSteps."
