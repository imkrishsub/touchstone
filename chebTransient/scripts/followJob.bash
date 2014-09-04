#!/bin/bash
qsub -W depend=afterany:$1.hopque01@hopque01 jobScriptUnfinishedCases.csh
