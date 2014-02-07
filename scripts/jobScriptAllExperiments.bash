#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=12:00:00
#PBS -N touchstone
#PBS -A m1041
#PBS -m ae

cd $PBS_O_WORKDIR

module load python

if [[ -z $p ]] ; then
  echo "please set p by running \'qsub -v p=x.xx ...\'"
  exit 1
fi

echo p=$p 
slope_ref=778.5
C_ref=7.624e6
lambda_ref=2.0
A_ref=4.6416e-25
codeFolder=../../code
CFL=400
dt=1e-6
maxXg=100.0
minXg=0.6936416184971098

mkdir -p vary_A
cd vary_A
for A in 4.6416e-24 2.1544e-24 1.0000e-24 4.6416e-25 2.1544e-25 1.0000e-25 \
  4.6416e-26 2.1544e-26 1.0000e-26
do
  prefix=A_$A
  commonArgs="--folder=. --p=$p --A=$A --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

mkdir -p vary_C
cd vary_C
C_CFL=50
for C in 3.812e8 7.624e7 7.624e6 7.624e5 3.812e5 2.541e5 1.906e5 \
  1.525e5 7.624e4 #3.812e4
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$C_CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
C_CFL=10
for C in 7.624e8 
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$C_CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

mkdir -p vary_slope
cd vary_slope
CFL=400
for slope in 500 1000 2500 5000 10000 50000
do
  prefix=slope_$slope
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
CFL=100
for slope in 100
do
  prefix=slope_$slope
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

mkdir -p vary_lambda_0
CFL=400
cd vary_lambda_0
for lambda_0 in 2.500e-1 5.000e-1 1.000e+0 2.000e+0 4.000e+0 8.000e+0 \
  1.600e+1 6.400e+1 1.280e+2
do
  prefix=lambda_0_$lambda_0
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
CFL=100
for lambda_0 in 3.125e-2 6.250e-2 
do
  prefix=lambda_0_$lambda_0
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  $codeFolder/scripts/runOneCase.bash "$commonArgs" $prefix $codeFolder &
done
cd ..

wait
