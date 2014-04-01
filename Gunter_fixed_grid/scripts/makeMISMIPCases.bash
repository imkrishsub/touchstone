#!/bin/bash

fileName=allMISMIPCases.txt
rm $fileName 2>/dev/null

slope_ref=778.5
C_ref=7.624e6
lambda_ref=2.0
dt=3e-4
xc=1.760


dxs=( 3.2 1.6 0.8 0.4 0.2 0.1 )
Nxs=( 551 1101 2201 4401 8801 17601 )
glpStrings=( "" "--useGLP" )
glpDirs=( nonGLP GLP)
commonArgs="--folder=. --linearSlope=$slope_ref --C=$C_ref --lambda_0=$lambda_ref --maxSteps=1000000"
# uncomment the following to include plotting
#commonArgs="$commonArgs --plot --plotContinuous" 

for ((resIndex=0; resIndex<=5; resIndex++))
do
  for p in 0.00 0.25 0.50 0.75 1.00
  do
    for ((glpIndex=0; glpIndex<=1; glpIndex++))
    do
      glpDir=${glpDirs[$glpIndex]}
      glpstring=${glpStrings[$glpIndex]}
      Nx=${Nxs[$resIndex]}
      dx=${dxs[$resIndex]}
      prevResult=none
      for A in 4.6416e-24 2.1544e-24 1.0000e-24 4.6416e-25 2.1544e-25 1.0000e-25 \
        4.6416e-26 2.1544e-26 1.0000e-26
      do
        prefix=A_${A}_adv
        args="$commonArgs --p=$p --A=$A --dt=$dt --deltaX=${dx}e-3 --Nx=$Nx --xc=$xc $glpString"
        echo $glpDir/p_$p/res_$dx >> $fileName
        echo $prefix >> $fileName
        echo $prevResult >> $fileName
        echo $args >> $fileName
        prevResult=${prefix}_final.pyda
      done
      for A in 2.1544e-26 4.6416e-26 1.0000e-25 2.1544e-25 4.6416e-25 1.0000e-24 \
        2.1544e-24 4.6416e-24
      do
        prefix=A_${A}_ret
        args="$commonArgs --p=$p --A=$A --dt=$dt --deltaX=${dx}e-3 --Nx=$Nx --xc=$xc $glpString"
        echo $glpDir/p_$p/res_$dx >> $fileName
        echo $prefix >> $fileName
        echo $prevResult >> $fileName
        echo $args >> $fileName
        prevResult=${prefix}_final.pyda
      done
    done
  done
done
