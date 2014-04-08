#!/bin/bash 

NSPES="01 02 03 04 05 06"
NRUNS="01 02 03 04 05 06 07 08 09 10"

if [ $# -gt 0 ] ; then
	NSPES="$1"
	NRUNS="$2"
fi

echo "NSPES: $NSPES"
echo "NRUNS: $NRUNS"

for n in $NSPES ; do
for i in $NRUNS ; do
	./fixedgrid $n 2>&1 | tee Output/run.$n.$i
done
done
