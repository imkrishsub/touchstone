#!/usr/bin/python
import os

caseFile = "allPolyMISMIPCases.txt"

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

for mismipIndex in range(70):
  for caseIndex in range(starts[mismipIndex],ends[mismipIndex]):
    AIndex = caseIndex-starts[mismipIndex]
    lineIndex=4*caseIndex
    dir = caseList[lineIndex]
    prefix = caseList[lineIndex+1]
    if(os.path.exists("%s/%s_final.pyda"%(dir,prefix))):
      print "%i %i: %s/%s final."%(mismipIndex, AIndex, dir,prefix)
      continue
    logFile="%s/%s_final.log"%(dir,prefix)
    if(os.path.exists(logFile)):
      error=False
      for line in open(logFile):
        if "blew" in line:
          print "%i %i: %s/%s failed in final."%(mismipIndex, AIndex, dir,prefix)
          error=True
          break
        if "Error" in line:
          print "%i %i: %s/%s failed in final."%(mismipIndex, AIndex, dir,prefix)
          error=True
          break
      if error:
        continue
  
    logFile="%s/%s_inProgress.log"%(dir,prefix)
    if(os.path.exists(logFile)):
      error=False
      for line in open(logFile):
        if "blew" in line:
          print "%i %i: %s/%s failed in inProgress."%(mismipIndex, AIndex, dir,prefix)
          error=True
          break
        if "Error" in line:
          print "%i %i: %s/%s failed in inProgress."%(mismipIndex, AIndex, dir,prefix)
          error=True
          break
      if not error:
        print "%i %i: %s/%s inProgress."%(mismipIndex, AIndex, dir,prefix)
