#!/bin/bash

find ready_for_analysis \( -name \*_1.fastq.gz -o -name \*[0-9][0-9].fastq.gz \
\) -print0 | xargs -0 dirname | parallel -j 6 -k -v ./second_queue_cleanup.sh
