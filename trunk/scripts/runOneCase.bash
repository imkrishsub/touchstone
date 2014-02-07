#!/bin/bash

commonArgs=$1 
prefix=$2
codeFolder=$3
inputFile=$4

echo running case: $prefix from $inputFile

#echo $commonArgs

initFile=${prefix}_init
looseFile=${prefix}_loose
strictFile=${prefix}_strict
finalFile=${prefix}_final
filePointer=${prefix}.pointer

commonCommand="aprun -n 1 python $codeFolder/main.py $commonArgs --filePointer=$filePointer"

restartFile=""
if [[ -f $filePointer ]] ; then
  restartFile=$(<$filePointer)
else
  if [[ -f $inputFile ]] ; then
    restartFile=$inputFile
  else
    echo warning: input file $inputFile not found.
  fi
fi

echo restarting from file: $restartFile

if [[ -z "$restartFile" ]] ; then
  echo init
  $commonCommand --outFile=$initFile.pyda \
    --eps_s=1e-3 --maxStep=0 --maxToleranceInner=1e-3 > $initFile.log
  if [ $? -ne 0 ]; then
    echo "init failed! Exiting."
    exit 1
  fi
  restartFile=$initFile.pyda
fi

if [ "$restartFile" != "$strictFile.pyda" ] && [ "$restartFile" != "$finalFile.pyda" ] ; then
  echo loose from $restartFile
  $commonCommand --outFile=$looseFile.pyda \
    --eps_s=1e-3 --maxStep=10000 --maxToleranceInner=1e-3 \
    --inFile=$restartFile --toleranceH=1e-1 --toleranceXg=1e-1 > $looseFile.log
  if [ $? -ne 0 ]; then
    echo "loose failed! Exiting."
    exit 1
  fi
 restartFile=$looseFile.pyda
fi

if [ "$restartFile" = "$looseFile.pyda" ] || [ "$restartFile" = "$strictFile.pyda" ] ; then
  echo strict from $restartFile
  $commonCommand --outFile=$strictFile.pyda \
    --eps_s=1e-8 --maxStep=10000 --maxToleranceInner=1e-5 \
    --inFile=$restartFile --toleranceH=1e-2 --toleranceXg=1e-2 \
    > $strictFile.log
  if [ $? -ne 0 ]; then
    echo "strict failed! Exiting."
    exit 1
  fi
 restartFile=$strictFile.pyda
fi

if [ "$restartFile" = "$strictFile.pyda" ] ; then
  echo final from $restartFile
  $commonCommand --outFile=$finalFile.pyda \
    --eps_s=1e-8 --maxStep=10000 --maxToleranceInner=1e-6 \
    --inFile=$restartFile --toleranceH=1e-2 --toleranceXg=1e-2 \
    > $finalFile.log
  if [ $? -ne 0 ]; then
    echo "final failed! Exiting."
    exit 1
  fi
fi

