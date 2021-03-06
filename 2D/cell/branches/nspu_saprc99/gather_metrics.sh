#!/bin/bash

MAX_SPUS=8

n=1
while [ $n -le $MAX_SPUS ] ; do

echo "
# Job name
#PBS -N fixedgrid
# Job queue
#PBS -q sdk21
# request 1 node
#PBS -l nodes=1:ppn=2
# Job wallclock time
#PBS -l walltime=00:15:00        

# Go to directory from which the job was submitted
cd \$PBS_O_WORKDIR

# Execute
./fixedgrid $n
" | qsub

((n++))

done
