#!/bin/bash

fileName=allCases.txt
rm $fileName

for p in 0.00 0.25 0.50 0.75 1.00
do

slope_ref=778.5
C_ref=7.624e6
lambda_ref=2.0
A_ref=4.6416e-25
dt=1e-6
maxXg=2.0
minXg=0.6936416184971098

C_init_adv = (1.994e6 7.624e6)
C_adv = (7.624e6 1.994e7)
CIndex = ${#C_adv[@]}

for ((index=0; i<${CImdex}; index++ )); 
do
  
  prefix=C_$C
  
  commonArgs="--folder=. --writeToSperateFile --fixedTimeStep --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref  --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  echo p_$p/vary_C >> $fileName
  echo $prefix >> $fileName
  echo $commonArgs >> $fileName
done



C_init_ret = (1.994e7 7.624e6)
C_ret = (7.624e6 1.994e6)

for C in 1.994e7 7.624e6 1.994e6
do
  prefix=C_$C
  commonArgs="--folder=. --p=$p --A=$A_ref --linearSlope=$slope_ref --C=$C --lambda_0=$lambda_ref --goalCFL=$CFL --initDt=$dt --maxXg=$maxXg --minXg=$minXg"
  echo p_$p/vary_C >> $fileName
  echo $prefix >> $fileName
  echo $commonArgs >> $fileName
done


