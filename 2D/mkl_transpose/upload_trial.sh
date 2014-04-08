#!/bin/bash

experiment="$1"
if [ -z "$experiment" ] ; then
  echo "Usage: $0 experiment_name"
  exit 1
fi

export PATH=$WORK/tau2/mic_linux/bin:$PATH

paraprof --pack $experiment.ppk
err=$?
if [ $err -eq 0 ] ; then 
  rm -f profile.* 
  rm -rf MULTI__*
  taudb_loadtrial -m "Experiment=$experiment:Application=fixedgrid" -n "fixedgrid-$experiment" -c stampede $experiment.ppk
fi

