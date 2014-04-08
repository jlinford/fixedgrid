#!/bin/bash

MAX_SPUS=6

n=1
while [ $n -le $MAX_SPUS ] ; do

ssh -f ps0$n "cd workspace/fixedgrid_3d/cell/trunk && ./dropbox_fixedgrid.sh $n"

((n++))

done
