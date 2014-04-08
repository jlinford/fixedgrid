#!/bin/bash

#for nprocs in 1 2 3 4 5 6 ; do
for nprocs in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ; do
	./fixedgrid $nprocs
	mv Output/METRICS_100.csv Output/METRICS_100_$nprocs.csv
done
