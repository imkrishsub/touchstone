#!/usr/bin/python

from optparse import OptionParser
import os
import subprocess 

parser = OptionParser()
options, args = parser.parse_args()

caseFile = args[0]
mismipIndex = int(args[1])

print 'caseFile=', caseFile 

caseList = [line.rstrip('\n') for line in open(caseFile)]

ACount = 17
for caseIndex in range(ACount):
  lineIndex=4*(mismipIndex+caseIndex)
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "case %s/%s already finished."%(dir,prefix)
    continue
  args = ["python", "code/scripts/runOneCase.py", caseFile, "%i"%(mismipIndex+caseIndex)]  
  
  print "running case %s/%s."%(dir,prefix)
  status = subprocess.call(args)
  if status != 0:
    print "runOneCase.py failed! Exiting."
    exit(status)  
  
