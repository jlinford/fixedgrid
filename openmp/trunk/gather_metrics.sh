#!/bin/bash

for nprocs in 1 2 3 4 5 6 7 ; do
	./fixedgrid $nprocs
	mv Output/METRICS_100.csv Output/METRICS_100_$nprocs.csv
done
