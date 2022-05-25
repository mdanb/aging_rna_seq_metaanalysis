#!/bin/bash
MAX_FILES=6
MAX_OCCUPIED_PERCENT=90
LOWER_THRESHOLD=30
DISK_SPACE=^/dev/nvme0n1p2

trap 'pkill -f first_queue.sh; pkill -f pigz; pkill -f second_queue.sh;         
pkill -f prefetch; pkill -f fastq-dump; pkill -f xargs; pkill -f kallisto;      
exit' INT

mkdir -p ready_for_analysis download_dir

if [ ! -f ./Caenorhabditis_elegans.WBcel235.cdna.all.index ]; then
	echo "Creating index file..."
	kallisto index -i Caenorhabditis_elegans.WBcel235.cdna.all.index \
	Caenorhabditis_elegans.WBcel235.cdna.all.fa
	echo "Done creating index file"
fi

# if you stop executing monitor_disk_space.sh, there might be still files that 
# haven't been mapped in the directory ready_for_analysis that should been 
# mapped. There will be also be large .gz files that will take up too much 
# space, and will cause disk space to overflow if you start downloading new 
# files. second_queue_cleanup.sh will make sure those remaining files are 
# mapped, and their corresponding .gz files deleted, so that you can continue
# downloading + mapping new files
./run_second_queue_cleanup.sh

# monitors the directory "ready_for_analysis" and maps *gz files that appear 
# there after first_queue.sh moves them there
./second_queue.sh &

while true; do
	percent_occupied=$(df -H | grep -E $DISK_SPACE | awk '{ print $5 }' \
	| sed 's/.$//')
	cat accession_list.txt | xargs -P$MAX_FILES -n3 ./first_queue.sh &
	while [[ $percent_occupied -lt $MAX_OCCUPIED_PERCENT ]]; do
		percent_occupied=$(df -H | grep -E $DISK_SPACE | awk \
		'{ print $5 }' | sed 's/.$//') 
		sleep 5
	done

	printf "surpassed $MAX_OCCUPIED_PERCENT %% disk space!\n"

	pkill -f first_queue.sh
	pkill -f fasterq-dump
	pkill -f prefetch
	find . -type d -name *tmp* | xargs rm -rf {} \;
	./first_queue_cleanup.sh &
	./run_second_queue_cleanup.sh &
	while [[ $percent_occupied -gt $LOWER_THRESHOLD ]]; do
		echo "$percent_occupied% currently occupied"
		percent_occupied=$(df -H | grep -E $DISK_SPACE | awk \
		'{ print $5 }' | sed 's/.$//')
		sleep 5
	done
done
