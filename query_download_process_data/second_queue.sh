#!/bin/bash

set -e

MONITORDIR="ready_for_analysis"
inotifywait --exclude _2 -m -r -e moved_to --format '%w%f' "${MONITORDIR}" | \
while read NEWFILE; do
	file_dir=$(dirname $NEWFILE)
	if [[ $NEWFILE =~ _1 ]]; then
		# wait until the other file for paired end is here (if paired end)
		while (( $(ls -1 $file_dir | wc -l) != 2 )); do
			sleep 1
		done
		{ set -e ; kallisto quant -i Caenorhabditis_elegans.WBcel235.cdna.all.index --plaintext -o $file_dir -t 12 $(ls $file_dir/*gz -1 | head -1) $(ls $file_dir/*gz -1 | tail -1) ; rm $file_dir/*gz ;  } &
	else
		{ set -e ; kallisto quant -i Caenorhabditis_elegans.WBcel235.cdna.all.index --plaintext -o $file_dir --single -t 12 -l 250 -s 30 $(ls $file_dir/*gz) ; rm $file_dir/*gz ;  } &
	fi
done


