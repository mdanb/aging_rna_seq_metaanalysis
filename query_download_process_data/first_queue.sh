#!/bin/bash

set -e # exit directly on pipeline error

if [ "$#" -ne 3 ]; then
    echo 'Number of arguments to first_queue must be 3!'
	echo 'Exiting with non 0 status'
	exit
fi   

srr=$1
bioproject=$2
library_layout=$3 

if [ -d ready_for_analysis/$bioproject/$srr ] ; then
	for file in $(find ready_for_analysis/$bioproject/$srr -type f); do
		if [[ $file =~ gz$ ]] || [[ $file =~ tsv ]]; then
			exit
		fi
	done
fi

# tmp file created from fasterq previously will no longer be used
# instead new one will be created
#if ! [ -d download_dir/${bioproject}_${library_layout}/$srr/*tmp* ] ; then

#if ! compgen -G "download_dir/${bioproject}_${library_layout}/$srr/*tmp*" > /dev/null; then
#	prefetch $srr -O download_dir/${bioproject}_${library_layout}/
#else
#	rm -rf download_dir/${bioproject}_${library_layout}/$srr/*tmp*
#fi


#if compgen -G "download_dir/${bioproject}_${library_layout}/$srr/*tmp*" > /dev/null; then
#	   rm -rf download_dir/${bioproject}_${library_layout}/$srr/*tmp*               
#fi

#if [ -f download_dir/${bioproject}_${library_layout}/$srr/*tmp$ ]; then
#	rm download_dir/${bioproject}_${library_layout}/$srr/*tmp$
#fi 

prefetch $srr -O download_dir/${bioproject}_${library_layout}/              

#cd download_dir/${bioproject}_${library_layout}/$srr
#temp_file=$(find . -maxdepth 1 -type f | egrep ".sra$")

#if (( $(find . -type f -maxdepth 1 | wc -l) == 1 )); then
#fasterq-dump $temp_file
#fi

cd download_dir/${bioproject}_${library_layout}/$srr
fasterq-dump ${srr}.sra
rm ${srr}.sra

#rm $temp_file

if [[ $library_layout == "PAIRED" ]]; then
	pigz --fast *_2.fastq &
	pigz --fast *_1.fastq &
else
	pigz --fast *.fastq &
fi

wait

mkdir -p ../../../ready_for_analysis/$bioproject/$srr
mv *gz ../../../ready_for_analysis/$bioproject/$srr
