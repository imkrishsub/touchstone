#!/usr/bin/python
import os

caseFile = "allMISMIPCases.txt"

caseList = [line.rstrip('\n') for line in open(caseFile)]

caseCount = 17*7*5*2
for caseIndex in range(caseCount):
  mismipIndex = caseIndex/17
  lineIndex=4*caseIndex
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "%i: %s/%s final."%(mismipIndex,dir,prefix)
    continue
  logFile="%s/%s_final.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i: %s/%s failed in final."%(mismipIndex,dir,prefix)
        error=True
        break
      if "Error" in line:
        error=True
        break
        print "%i: %s/%s failed in final."%(mismipIndex,dir,prefix)
    if error:
      continue

  logFile="%s/%s_inProgress.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i: %s/%s failed in inProgress."%(mismipIndex,dir,prefix)
        error=True
        break
      if "Error" in line:
        error=True
        break
        print "%i: %s/%s failed in inProgress."%(mismipIndex,dir,prefix)
    if error:
      continue
    print "%i: %s/%s inProgress."%(mismipIndex,dir,prefix)
