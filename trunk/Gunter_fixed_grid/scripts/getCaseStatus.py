#!/usr/bin/python
import os

caseFile = "allMISMIPCases.txt"

caseList = [line.rstrip('\n') for line in open(caseFile)]

caseCount = 14*7*5*2
for caseIndex in range(caseCount):
  mismipIndex = caseIndex/14
  CIndex = caseIndex-14*mismipIndex
  lineIndex=4*caseIndex
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "%i %i: %s/%s final."%(mismipIndex, CIndex, dir,prefix)
    continue
  logFile="%s/%s_final.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i %i: %s/%s failed in final."%(mismipIndex, CIndex, dir,prefix)
        error=True
        break
      if "Error" in line:
        print "%i %i: %s/%s failed in final."%(mismipIndex, CIndex, dir,prefix)
        error=True
        break
    if error:
      continue

  logFile="%s/%s_inProgress.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i %i: %s/%s failed in inProgress."%(mismipIndex, CIndex, dir,prefix)
        error=True
        break
      if "Error" in line:
        print "%i %i: %s/%s failed in inProgress."%(mismipIndex, CIndex, dir,prefix)
        error=True
        break
    if not error:
      print "%i %i: %s/%s inProgress."%(mismipIndex, CIndex, dir,prefix)
