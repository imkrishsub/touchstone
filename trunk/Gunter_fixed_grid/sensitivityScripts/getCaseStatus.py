#!/usr/bin/python
import os

caseFile = "allSensitivityCases.txt"

caseList = [line.rstrip('\n') for line in open(caseFile)]

caseCount = len(caseList)/4

exptFirstCase = []
for caseIndex in range(caseCount):
  lineIndex=4*caseIndex

exptFirstCase.append(caseCount)

exptIndex = 0
caseNumber = 0
for caseIndex in range(caseCount):
  lineIndex=4*caseIndex
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  prevResult=caseList[lineIndex+2]
  if(prevResult == "none"):
    caseNumber = 0
  if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
    print "%i %i: %s/%s final."%(exptIndex, caseNumber, dir,prefix)
    continue
  logFile="%s/%s_final.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i %i: %s/%s failed in final."%(exptIndex, caseNumber, dir,prefix)
        error=True
        break
      if "Error" in line:
        print "%i %i: %s/%s failed in final."%(exptIndex, caseNumber, dir,prefix)
        error=True
        break
    if error:
      continue

  logFile="%s/%s_inProgress.log"%(dir,prefix)
  if(os.path.exists(logFile)):
    error=False
    for line in open(logFile):
      if "blew" in line:
        print "%i %i: %s/%s failed in inProgress."%(exptIndex, caseNumber, dir,prefix)
        error=True
        break
      if "Error" in line:
        print "%i %i: %s/%s failed in inProgress."%(exptIndex, caseNumber, dir,prefix)
        error=True
        break
    if not error:
      print "%i %i: %s/%s inProgress."%(exptIndex, caseNumber, dir,prefix)

  if(prevResult == "none"):
    exptIndex += 1
  caseNumber += 1

