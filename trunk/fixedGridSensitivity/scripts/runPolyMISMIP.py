#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 
import numpy

parser = OptionParser()
options, args = parser.parse_args()

caseFile = args[0]
mismipIndex = int(args[1])

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]
exptCount = int(caseList[0])
exptIndex = 0
lineIndex = 1
startCase = 0
ACount = int(caseList[lineIndex])
while(exptIndex < mismipIndex):
  lineIndex += 4*ACount+1
  startCase += ACount
  exptIndex += 1
  ACount = int(caseList[lineIndex])

endCase = startCase+ACount

print exptCount, exptIndex, startCase, endCase, lineIndex

lineIndex += 1 # skip the A count
for caseIndex in range(startCase,endCase):
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "case %s/%s already finished."%(dir,prefix)
    lineIndex += 4
    continue
  args = ["python", "code/scripts/runOneCase.py", caseFile, "%i"%(caseIndex), "%i"%(lineIndex)]  
  
  print "running case %s/%s."%(dir,prefix)
  status = subprocess.call(args)
  if status != 0:
    print "runOneCase.py failed! Exiting."
    exit(status)  
  lineIndex += 4
  
