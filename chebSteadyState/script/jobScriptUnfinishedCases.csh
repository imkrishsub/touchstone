#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=48
#PBS -l walltime=24:00:00
#PBS -N touchstone
#PBS -A m1041
#PBS -m ae

module load taskfarmer
module load python

cd $PBS_O_WORKDIR

set caseFile = unfinishedCases.txt
set logDir = logs$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

set caseLines = `wc -l < $caseFile`
echo caseLines: $caseLines
@ caseCount = ( $caseLines / 3 )

echo caseCount: $caseCount

set script = "python code/runOneCase.py"
echo script: $script

tf -t $caseCount -n 2 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID
