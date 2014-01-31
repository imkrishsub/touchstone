#!/bin/bash
for p in 0.00 0.25 0.50 0.75 1.00
do
  cd p_$p
  qsub -v p=$p ../code/scripts/jobScriptAllExperiments.bash
  cd ..
done
