#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=48
#PBS -l walltime=2:00:00
#PBS -N touchstone
#PBS -A m1041
#PBS -m ae
#PBS -V

cd $PBS_O_WORKDIR
#exptFile=experiments1.txt

module load python
module load taskfarmer

if [[ -z $caseFile ]] ; then
  echo "A value for expFile must be set by running \'qsub -v caseFile=xxx.txt ...\' Exiting."
  exit 1
fi

echo caseFile=$caseFile 

IFS=$'\r\n' caseList=($(cat $caseFile))

caseLines=${#caseList[@]}
echo caseLines: $caseLines
caseCount=$(( caseLines / 3 ))

echo caseCount: $caseCount

tf -t $caseCount -o touchstoneCase.%t ./code/scripts/runOneCase.bash $caseFile \$TF_TASKID
