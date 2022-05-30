#!/bin/bash
#for file in $(find download_dir | grep \.sra$); do
#	rm $file	
#done

find download_dir/ -name *.sra -exec rm -rf {} \;

for file in $(find download_dir -name *.gz); do
	bioproject=$(dirname $file | awk -F / '{ print $2' } | awk -F _ '{ print
	$1 }')
	srr=$(dirname $file | awk -F / '{ print $2' }) 
	#bioproject=$(echo $file | cut -d / -f 2 | cut -d _ -f 1)
	#srr=$(echo $file | cut -d / -f 3)
    mkdir -p ready_for_analysis/$bioproject/$srr	
	mv $file ready_for_analysis/$bioproject/$srr
done

for file in $(find download_dir -name *fastq*); do
	bioproject=$(dirname $file | awk -F / '{ print $2' } | awk -F _ '{ print    
    $1 }')
	srr=$(dirname $file | awk -F / '{ print $2' })
	#bioproject=$(echo $file | cut -d / -f 2 | cut -d _ -f 1)
	#srr=$(echo $file | cut -d / -f 3)
	if [[ $file =~ PAIRED ]]; then
		if [[ $file =~ _2 ]] || [[ $file =~ _1 ]]; then
			{ set -e; pigz --fast $file ; mkdir -p 
			ready_for_analysis/$bioproject/$srr ; mv $file.gz ready_for_analysis/$bioproject/$srr ; }&
		else
			echo "$file" >> trouble_files
			rm $file
			sed -i "/$(echo $file | rev | cut -d / -f1 | rev)/d" accession_list
		fi
	else
		{ set -e; pigz --fast $file ; mkdir -p ready_for_analysis/$bioproject/$srr ; mv $file.gz ready_for_analysis/$bioproject/$srr ; }&
	fi
done
