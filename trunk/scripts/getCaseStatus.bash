#!/bin/bash

shopt -s nullglob
#initFiles=p*/vary_*/*init.log
initFiles=p*/*/*init.log
#echo $initFiles
for initFile in $initFiles
do
  IFS='/' read pFolder vFolder fileName <<< "$initFile"
  IFS='_' read var value ending <<< "$fileName"
  expt=${var}_${value}
  prefix=$pFolder/$vFolder/$expt
  if [[ -f ${prefix}_final.log ]] ; then
    echo $prefix: final
  elif [[ -f ${prefix}_strict.log ]] ; then
    echo $prefix: strict
  elif [[ -f ${prefix}_loose.log ]] ; then
    echo $prefix: loose
  else
    echo $prefix: init
  fi

done
