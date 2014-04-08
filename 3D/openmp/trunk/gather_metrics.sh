#!/bin/bash

for n in 1 2 3 4 5 6 7 8 ; do
	echo -n "Running fixedgrid $n..."
	./fixedgrid $n 2>&1 > Output/fixedgrid_$n.out
	echo " done!"
done
