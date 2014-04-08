#!/bin/bash

NUMBER="all"

if [ "${NUMBER}${1}" != "${NUMBER}" ] ; then
	NUMBER="$1"
fi

PREFIX="OUT_solution_O3_100_00"
SUFFIX=".000"
REF_DIR="reference_output"

echo "Verifying against $REF_DIR..."

if [ "$NUMBER" == "all" ] ; then
	for file in `ls Output` ; do
		f=`echo $file | cut -d/ -f2`
		diff Output/$f $REF_DIR/$f
	done
else
	f="${PREFIX}${NUMBER}${SUFFIX}"
	diff Output/$f $REF_DIR/$f
fi

echo "Done."
