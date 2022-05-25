#!/bin/bash

set -e
if [ $# -ne 1 ] ; then
    echo "Wrong number of arguments"
    exit 1
else
    echo "I got $1 as my argument"
fi

file_dir=$1
count=$(ls $file_dir/*.gz | wc -l)
if [[ $count -eq 2 ]]; then
	kallisto quant -i Caenorhabditis_elegans.WBcel235.cdna.all.index \
	--plaintext -o $file_dir -t 12 $(ls $file_dir/*gz | head -1) \
	$(ls $file_dir/*gz | tail -1)
	rm $file_dir/*gz # paired end
elif [[ $count -eq 1 ]]; then
	kallisto quant -i Caenorhabditis_elegans.WBcel235.cdna.all.index \
	--plaintext -o $file_dir --single -t 12 -l 250 -s 30 $(ls $file_dir/*gz)
	rm $file_dir/*gz # single end
else
	echo "More than two files were found when running kallisto"
	echo "Exiting with non zero status"
	exit 1
fi
