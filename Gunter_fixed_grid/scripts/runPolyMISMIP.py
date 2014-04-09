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

ACount_p = [37, 37, 37, 67, 67]
starts = []
ends = []
start = 0
for resIndex in range(7):
  for pIndex in range(5):
    for glpIndex in range(2):
      starts += [start]
      start += ACount_p[pIndex]
      ends += [start]

for caseIndex in range(starts[mismipIndex], ends[mismipIndex]):
  lineIndex=4*caseIndex
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "case %s/%s already finished."%(dir,prefix)
    continue
  args = ["python", "code/scripts/runOneCase.py", caseFile, "%i"%(caseIndex)]  
  
  print "running case %s/%s."%(dir,prefix)
  status = subprocess.call(args)
  if status != 0:
    print "runOneCase.py failed! Exiting."
    exit(status)  
  
