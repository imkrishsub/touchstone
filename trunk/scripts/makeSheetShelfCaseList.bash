#!/bin/bash

fileName=allCases.txt
rm $fileName

for p in 0.00 0.25 0.50 0.75 1.00
do

slope_ref=778.5
C_ref=7.624e6
lambda_ref=2.0
A_ref=2.1544e-25
dt=1e-6
maxXg=2.112
minXg=0.6936416184971098

#CFL=200
#for A in 4.6416e-24 2.1544e-24 1.0000e-24 4.6416e-25 2.1544e-25 1.0000e-25 \
#  4.6416e-26 2.1544e-26 1.0000e-26
#do
#  prefix=A_$A
#  commonArgs="--folder=. --p=$p --A=$A --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
###  echo p_$p/vary_A >> $fileName
#  echo $prefix >> $fileName
#  echo $commonArgs >> $fileName
#done

CFL=50
for C in 3.160e5 5.008e5 7.938e5 1.258e6 1.994e6 3.160e6 5.008e6 7.624e6 1.258e7 1.994e7 3.160e7
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  echo p_$p/vary_C >> $fileName
  echo $prefix >> $fileName
  echo $commonArgs >> $fileName
done

CFL=200
for slope in 155.70e+2 778.50e+1 389.25e+1 311.40e+1 233.55e+1 155.70e+1 778.50e+0 584.00e+0
do
  prefix=slope_$slope
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  echo p_$p/vary_slope >> $fileName
  echo $prefix >> $fileName
  echo $commonArgs >> $fileName
done

#CFL=400
#for lambda_0 in 3.125e-2 6.250e-2 1.250e-1 2.500e-1 5.000e-1 1.000e+0 2.000e+0 4.000e+0 8.000e+0 1.600e+1 3.200e+1 6.400e+1 1.280e+2
#do
#  prefix=lambda_0_$lambda_0
#  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
#  echo p_$p/vary_lambda_0 >> $fileName
#  echo $prefix >> $fileName
#  echo $commonArgs >> $fileName
#done

#CFL=100
#for lambda_0 in 3.125e-2 6.250e-2 
#do
#  prefix=lambda_0_$lambda_0
#  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_0 --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
#  echo p_$p/vary_lambda_0 >> $fileName
#  echo $prefix >> $fileName
#  echo $commonArgs >> $fileName
#done

#done

#for p in 0.00 0.25 0.50 0.75 1.00
#do
CFL=25
for W in 0.5 1 2 3 4 5 10 20 100
do
  prefix=W_$W
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg --useChannelWidth --W0=$W --Wx=0.0"
  echo p_$p/vary_W >> $fileName
  echo $prefix >> $fileName
  echo $commonArgs >> $fileName
done
done
