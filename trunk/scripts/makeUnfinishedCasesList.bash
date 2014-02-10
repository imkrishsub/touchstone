#!/bin/bash

inCaseFile=allCases.txt
outCaseFile=unfinishedCases.txt
IFS=$'\r\n' caseList=($(cat $inCaseFile))

caseLines=${#caseList[@]}
echo caseLines: $caseLines
caseCount=$(( caseLines / 3 ))

rm $outCaseFile

lineIndex=0
caseIndex=0
for index in $(seq 1 $caseCount)
do
  dir="${caseList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  prefix="${caseList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  commonArgs="${caseList[$lineIndex]}"
  let lineIndex=$lineIndex+1

  logPrefix=$dir/$prefix
  if [[ -f ${logPrefix}_final.pyda ]] ; then
    echo $index skipping finished case $logPrefix
  else
    echo $index adding unfinished case $caseIndex to $logPrefix
    echo $dir >> $outCaseFile
    echo $prefix >> $outCaseFile
    echo $commonArgs >> $outCaseFile
    let caseIndex=$caseIndex+1
  fi

done
