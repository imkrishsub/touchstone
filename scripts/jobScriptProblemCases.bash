#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=02:00:00
#PBS -N touchstoneTroubleCases
#PBS -A m1041
#PBS -m ae

cd $PBS_O_WORKDIR

slope_ref=778.5
C_ref=7.624e6
lambda_ref=2.0
A_ref=4.6416e-25
codeFolder=../../code
CFL=200
dt=1e-6
maxXg=100.0
minXg=0.6936416184971098

p=0.50
cd p_$p

CFL=5
cd vary_C
for C in 7.624e8
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

CFL=100
cd vary_slope
for slope in 100
do
  prefix=slope_$slope
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

cd ..
p=0.75
cd p_$p

CFL=100
cd vary_lambda_0
for lambda_0 in 3.125e-2 6.250e-2 
do
  prefix=lambda_0$lambda_0
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

cd ..
p=1.00
cd p_$p

CFL=100
cd vary_A
for A in 4.6416e-24 
do
  prefix=A_$A
  commonArgs="--folder=. --p=$p --A=$A --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

CFL=5
cd vary_C
for C in 7.624e8
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

CFL=100
cd vary_slope
for slope in 100
do
  prefix=slope_$slope
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

CFL=100
cd vary_lambda_0
for lambda_0 in 3.125e-2 6.250e-2
do
  prefix=lambda_0$lambda_0
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

wait
