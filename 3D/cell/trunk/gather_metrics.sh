#!/bin/bash

MAX_SPUS=16

n=0
while [ $n -le $MAX_SPUS ] ; do

if [ $n -lt 10 ] ; then
	zeros="00"
else
	zeros="0"
fi

echo "
# Job name
#PBS -N fixedgrid_$n
# Job queue
#PBS -q serial
# request 1 node
#PBS -l nodes=node$zeros$n
# Job wallclock time
#PBS -l walltime=8:0:0
#PBS -l cput=8:0:0

# Go to directory from which the job was submitted
cd \$PBS_O_WORKDIR

# Execute
./fixedgrid $n
" | qsub

((n++))

done
