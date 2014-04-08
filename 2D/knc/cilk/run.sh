#!/bin/bash

MIC=$MIC_ENV_PREFIX

export ${MIC}_TAU_VERBOSE=1

export ${MIC}_KMP_AFFINITY=compact:verbose

./fixedgrid

