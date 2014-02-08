#PBS -S /bin/csh
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=00:05:00
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

set script = "python runOneCase.py"
echo script: $script

tf -t $caseCount -n 1 -o $logDir/touchstoneCase%t.out -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID
