#!/bin/bash
setNumber=1
fileName=experiments$setNumber.txt
while [[ -f $fileName ]]
do
  qsub -v exptFile=$fileName code/scripts/jobScriptExperimentsFromList.bash
  let setNumber=$setNumber+1
done
