#!/bin/bash
slope=778.5
p=$1
C=7.624e6
lambda_0=2.0
codeFolder=../../code
dt=1e-6
CFL=$2
maxXg=100.0

prevResult=""
for A in 4.6416e-24 2.1544e-24 1.0000e-24 4.6416e-25 2.1544e-25 1.0000e-25 4.6416e-26 2.1544e-26 1.0000e-26
do
  prefix=A_$A
  commonArgs="--folder=. --p=$p --A=$A --linearSlope=$slope --C=$C --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --plot --plotContinuous"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder $prevResult
  if [ $? -ne 0 ]; then
    echo "runOneCase.bash failed! Exiting."
    exit 1
  fi

  prevResult=${prefix}_final.pyda
done

