#!/bin/bash

exptFile=experiments.txt
IFS=$'\r\n' exptList=($(cat $exptFile))

np=24

exptLines=${#exptList[@]}
echo exptLines: $exptLines
maxExpt=$(( exptLines / 3 ))

rm experiments?.txt

lineIndex=0
exptFileNumber=1
exptNumber=1
for index in $(seq 1 $maxExpt)
do
  echo $index $exptFileNumber $exptNumber
  dir="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  prefix="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  commonArgs="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1

  logPrefix=$dir/$prefix
  if [[ -f ${logPrefix}_final.pyda ]] ; then
    echo skipping finished case: $logPrefix
  else
    fileName=experiments$exptFileNumber.txt
    echo $dir >> $fileName
    echo $prefix >> $fileName
    echo $commonArgs >> $fileName

    let exptNumber=$exptNumber+1
    if (($exptNumber > $np)) ; then
      let exptFileNumber=$exptFileNumber+1
      exptNumber=1
    fi
  fi

done
