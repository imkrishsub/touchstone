#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=3:20:00
#PBS -N FixedGridMISMIP
#PBS -A m1041
#PBS -m ae

module load taskfarmer
module load python

cd $PBS_O_WORKDIR
set queue=edique02
echo $PBS_JOBID@$queue

set id=`qsub -W depend=afterany:$PBS_JOBID@$queue jobScriptMISMIP.csh`
echo submitted the next job: $id

mkdir -p logs
mv FixedGrid* logs

set finishedCount = `./getMISMIPStatus.bash | wc -l`
echo $finishedCount finished experiments will be copied to the results folder.
./copyResults.bash
@ remainingCount = ( 70 - $finishedCount )
echo $remainingCount unfinished experiments are still in progress.
set failedCount = `./getCaseStatus.bash | grep failed | wc -l`
echo $failedCount unfinished experiments have failed and need attention.

set caseFile = allMISMIPCases.txt
set logDir = logs/$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

set caseLines = `wc -l < $caseFile`
echo caseLines: $caseLines
set CCount = 14
@ linesPerExpt = ( 4 * $CCount )
@ caseCount = ( $caseLines / $linesPerExpt )

echo caseCount: $caseCount

set script = "python code/scripts/runMISMIP.py"
echo script: $script

tf -t $caseCount -n 1 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID
