#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=2:00:00
#PBS -N touchstone
#PBS -A m1041
#PBS -m ae

cd $PBS_O_WORKDIR
#exptFile=experiments1.txt

module load python
if [[ -z $exptFile ]] ; then
  echo "A value for expFile must be set by running \'qsub -v exptFile=xxx.txt ...\' Exiting."
  exit 1
fi

echo exptFile=$exptFile 

np=24
codeFolder=../../code

IFS=$'\r\n' exptList=($(cat $exptFile))

exptLines=${#exptList[@]}
echo exptLines: $exptLines
maxExpt=$(( exptLines / 3 ))
echo maxExpt: $maxExpt
lastExpt=$np
if (($lastExpt > $maxExpt)) ; then
  lastExpt=$maxExpt
fi

echo lastExpt: $lastExpt

lineIndex=0
for index in $(seq 1 $lastExpt)
do
  echo $index
  dir="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  mkdir -p $dir
  cd $dir
  prefix="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  commonArgs="${exptList[$lineIndex]}"
  let lineIndex=$lineIndex+1
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
  #echo $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder
  cd ../..
done

wait
