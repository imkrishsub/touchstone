#!/bin/bash

shopt -s nullglob
initFiles=p*/*/*init.log
for initFile in $initFiles
do
  IFS='/' read pFolder vFolder fileName <<< "$initFile"
  IFS='_' read var value ending <<< "$fileName"
  if [[ "$var" = "lambda" ]] ; then 
    IFS='_' read var extra value ending <<< "$fileName"
    var=${var}_$extra
  fi
  expt=${var}_${value}
  prefix=$pFolder/$vFolder/$expt
  if [[ -f ${prefix}_final.pyda ]] ; then
    mv ${prefix}_final.pyda ${prefix}_loose.pyda
    echo ${expt}_loose.pyda > ${prefix}.pointer
    rm ${prefix}_final.log ${prefix}_strict.*
    echo $prefix: final
  elif [[ -f ${prefix}_strict.log ]] ; then
    mv ${prefix}_strict.pyda ${prefix}_loose.pyda
    echo ${expt}_loose.pyda > ${prefix}.pointer
    rm ${prefix}_final.log ${prefix}_strict.log
    echo $prefix: strict
  fi
done
