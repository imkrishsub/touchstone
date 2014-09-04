#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 

parser = OptionParser()
options, args = parser.parse_args()

caseFile = args[0]
exptIndex = int(args[1])

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]

caseCount = len(caseList)/4

exptFirstCase = []
for caseIndex in range(caseCount):
  lineIndex=4*caseIndex
  prevResult=caseList[lineIndex+2]
  if(prevResult == "none"):
    exptFirstCase.append(caseIndex)

exptFirstCase.append(caseCount)

for caseIndex in range(exptFirstCase[exptIndex],exptFirstCase[exptIndex+1]):
  lineIndex=4*caseIndex
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "case %s/%s already finished."%(dir,prefix)
    continue
  args = ["python", "code/scripts/runOneCase.py", caseFile, "%i"%(caseIndex), "%i"%(lineIndex)]  
  
  print "running case %s/%s."%(dir,prefix)
  status = subprocess.call(args)
  if status != 0:
    print "runOneCase.py failed! Exiting."
    exit(status)  
  

