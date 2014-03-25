#PBS -S /bin/csh
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=24:00:00
#PBS -N FixedGridMISMIP
#PBS -A m1041
#PBS -m ae

module load taskfarmer
module load python

cd $PBS_O_WORKDIR

set caseFile = allMISMIPCases.txt
set logDir = logs$PBS_JOBID
echo logDir:$logDir
mkdir -p $logDir
echo caseFile: $caseFile

set caseLines = `wc -l < $caseFile`
echo caseLines: $caseLines
ACount=17
let linesPerExpt=4*$ACount
@ caseCount = ( $caseLines / linesPerExpt )

echo caseCount: $caseCount

set script = "python code/scripts/runMISMIP.py"
echo script: $script

tf -t $caseCount -n 1 -o $logDir/touchstoneCase%t.out \
  -e $logDir/touchstoneCase%t.err $script $caseFile \$TF_TASKID
