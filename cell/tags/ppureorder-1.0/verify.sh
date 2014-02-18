#!/bin/bash

NUMBER="all"

if [ "${NUMBER}${1}" != "${NUMBER}" ] ; then
	NUMBER="$1"
fi

VERIFY_CMD=./compare

PREFIX="OUT_solution_O3_08_000"
SUFFIX=".000"
REF_DIR="ref_out"

echo "Verifying against $REF_DIR..."

if [ "$NUMBER" == "all" ] ; then
	for file in `ls Output` ; do
		f=`echo $file | cut -d/ -f2`
		if [ "$f" != "METRICS_100.csv" ] ; then
			echo "$f"
			$VERIFY_CMD Output/$f $REF_DIR/$f
		fi
	done
else
	f="${PREFIX}${NUMBER}${SUFFIX}"
	$VERIFY_CMD Output/$f $REF_DIR/$f
fi

echo "Done."
