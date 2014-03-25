#!/bin/bash

shopt -s nullglob
initFiles=p*/*/*init.log
for initFile in $initFiles
do
  IFS='/' read pFolder resFolder fileName <<< "$initFile"
  IFS='_' read var value advRet ending <<< "$fileName"
  if [[ "$var" = "lambda" ]] ; then 
    IFS='_' read var extra value advRet ending <<< "$fileName"
    var=${var}_$extra
  fi
  expt=${var}_${value}_${advRet}
  prefix=$pFolder/$resFolder/$expt
  if [[ -f ${prefix}_final.log ]] ; then
    logFile=${prefix}_final.log
    echo $prefix: final
  elif [[ -f ${prefix}_strict.log ]] ; then
    logFile=${prefix}_strict.log
    echo $prefix: strict
  elif [[ -f ${prefix}_loose.log ]] ; then
    logFile=${prefix}_loose.log
    echo $prefix: loose
  else
    logFile=${prefix}_init.log
    echo $prefix: init
  fi
  result=`grep blew "$logFile"`
  if [[ -n $result ]] ; then
    echo $prefix failed! Try reducing CFL.
  fi
  result=`grep Error "$logFile"`
  if [[ -n $result ]] ; then
    echo $prefix failed! Try reducing CFL.
  fi

done
