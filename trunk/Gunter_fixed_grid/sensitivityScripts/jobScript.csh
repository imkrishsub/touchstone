#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=3:20:00
#PBS -N FixedGridSensitivity
#PBS -A m1041
#PBS -m ae

module load taskfarmer
module load python

cd $PBS_O_WORKDIR
set queue=edique02
echo $PBS_JOBID@$queue

set id=`qsub -W depend=afterany:$PBS_JOBID@$queue jobScript.csh`
echo submitted the next job: $id

mkdir -p logs
mv FixedGrid* logs

set finishedCount = `./code/sensitivityScripts/getCaseStatus.bash | wc -l`
echo $finishedCount finished experiments will be copied to the results folder.
./copyResults.bash
@ remainingCount = (280 - $finishedCount )
echo $remainingCount unfinished experiments are still in progress.
set failedCount = `./getCaseStatus.bash | grep failed | wc -l`
echo $failedCount unfinished experiments have failed and need attention.

set caseFile = allSensitivityCases.txt
set logDir = logs/$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

exptCount=280

echo exptCount: $exptCount

set script = "python code/scripts/runSensitivityExpt.py"
echo script: $script

tf -t $exptCount -n 1 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID

