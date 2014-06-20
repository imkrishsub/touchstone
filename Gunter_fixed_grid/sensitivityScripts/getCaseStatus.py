#!/usr/bin/python
import os

caseFile = "allSensitivityCases.txt"

caseList = [line.rstrip('\n') for line in open(caseFile)]

caseCount = len(caseList)/4


exptFirstCase = []
for caseIndex in range(caseCount):
  lineIndex=4*caseIndex
  prevResult=caseList[lineIndex+2]
  if(prevResult == "none"):
    exptFirstCase.append(caseIndex)

exptFirstCase.append(caseCount)

exptCount = len(exptFirstCase)-1
for exptIndex in range(exptCount):
  caseCount = exptFirstCase[exptIndex+1] - exptFirstCase[exptIndex]
  lastCase = -1
  for caseIndex in range(caseCount):
    lineIndex=4*(caseIndex + exptFirstCase[exptIndex])
    dir = caseList[lineIndex]
    prefix = caseList[lineIndex+1]
    exists = False
    if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
      stage = 'final'
      lastCase = caseIndex
    elif(os.path.exists("%s/%s_inProgress.pyda"%(dir,prefix))):
      stage = 'inProgress'
      lastCase = caseIndex
    else:
      break

  if lastCase == -1:
    lineIndex=4*(exptFirstCase[exptIndex])
    dir = caseList[lineIndex]
    prefix = caseList[lineIndex+1]
    print "expt %i: %s/%s not started."%(exptIndex,dir,prefix)
    continue

  lineIndex=4*(lastCase + exptFirstCase[exptIndex])
  dir = caseList[lineIndex]
  prefix = caseList[lineIndex+1]
  logFile="%s/%s_%s.log"%(dir,prefix,stage)
  if(not os.path.exists(logFile)):
    print "expt %i, case %i: %s/%s log file not found."%(exptIndex,lastCase,
      dir,prefix)
    continue
  error=False
  for line in open(logFile):
    if "blew" in line:
      error=True
      break
    if "Error" in line:
      error=True
      break
  if error:
    print "expt %i, case %i: %s/%s failed in %s."%(exptIndex,lastCase,
      dir,prefix,stage)
    continue
  print "expt %i, case %i: %s/%s %s."%(exptIndex,lastCase,
    dir,prefix,stage)

