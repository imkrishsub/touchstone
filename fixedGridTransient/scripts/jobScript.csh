#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=144
#PBS -l walltime=30:00
#PBS -N FixedGridTransient
#PBS -A m1041
#PBS -m ae

module load taskfarmer
module load python

cd $PBS_O_WORKDIR
set queue=edique02
echo $PBS_JOBID@$queue

#set id=`qsub -W depend=afterany:$PBS_JOBID@$queue jobScript.csh`
#echo submitted the next job: $id

mkdir -p logs
mv FixedGrid* logs

set exptCount = 140

echo exptCount: $exptCount

set caseFile = allTransientCases.txt
set logDir = logs/$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

set script = "python code/scripts/runOneCase.py"
echo script: $script

tf -t $exptCount -n 6 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID

