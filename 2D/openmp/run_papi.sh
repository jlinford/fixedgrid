#!/bin/bash

MIC=$MIC_ENV_PREFIX

for metric in `$TACC_PAPI_INC/../bin/papi_avail | grep Yes | awk '{print $1}'`
do
  export ${MIC}_TAU_METRICS=TIME:$metric
  ./run.sh
done

