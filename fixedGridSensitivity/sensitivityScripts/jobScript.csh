#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=288
#PBS -l walltime=30:00
#PBS -N FixedGridSensitivity
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


set finishedCount = `aprun -n 1 python code/sensitivityScripts/getCaseStatus.py | wc -l`
# subtract extra line because of aprun output
@ finishedCount = ( $finishedCount - 1 )
echo $finishedCount finished experiments.
@ remainingCount = ( 280 - $finishedCount )
echo $remainingCount unfinished experiments are still in progress.
set failedCount = `aprun -n 1 python code/sensitivityScripts/getCaseStatus.py | grep failed | wc -l`
# subtract extra line because of aprun output
@ failedCount = ( $failedCount - 1 )
echo $failedCount unfinished experiments have failed and need attention.

set caseFile = allSensitivityCases.txt
set logDir = logs/$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

set exptCount = 280

echo exptCount: $exptCount

set script = "python code/sensitivityScripts/runSensitivityExpt.py"
echo script: $script

tf -t $exptCount -n 12 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID

