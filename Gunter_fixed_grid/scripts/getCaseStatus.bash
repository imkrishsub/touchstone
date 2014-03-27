#!/bin/bash

shopt -s nullglob
initFiles=p*/*/*inProgress.log
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
  else
    logFile=${prefix}_inProgress.log
    echo $prefix: in progress
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
