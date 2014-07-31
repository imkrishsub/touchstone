#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 

parser = OptionParser()
options, args = parser.parse_args()

caseFile = args[0]
caseIndex = int(args[1])

lineIndex = 4*caseIndex
#lineIndex = int(args[2])

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]

dir=caseList[lineIndex]
prefix=caseList[lineIndex+1]
inputFile=caseList[lineIndex+2]
commonArgs=caseList[lineIndex+3].split()

print "running case: %s/%s"%(dir,prefix)
print "reading input from", inputFile

if(not os.path.exists(dir)):
  os.makedirs(dir)

oldPath = os.getcwd()

os.chdir(dir)

filePointer="%s.pointer"%prefix

commonArgs = ["python", "%s/code/Solver.py"%oldPath] \
  + commonArgs \
  + ["--folder=.", "--filePointer=%s"%filePointer]

restartFile="none"
if os.path.exists(filePointer):
  lines = [line.rstrip('\n') for line in open(filePointer)]
  if len(lines) > 0 and os.path.exists(lines[0]) and (os.stat(lines[0]).st_size > 0):
    restartFile=lines[0]
    print "restating from", restartFile
if(restartFile == "none"):
  if os.path.exists(inputFile):
    print "stating from input file", inputFile
    restartFile=inputFile
  else:
    print "exiting: input file", inputFile, "not found."
    exit(1)

logFile = open("%s.log"%prefix,'w')
errFile = open("%s.err"%prefix,'w')
args = commonArgs + ["--outFile=%s.pyda"%prefix, "--inFile=%s"%restartFile]
status = subprocess.call(args, stdout=logFile, stderr=errFile)
logFile.close()
errFile.close()
if status == 40:
  print "failed to read input file. Try modifying the file pointer to restart from an earlier time step."

if status != 0:
  print "run failed! Exiting."
  exit(status)
