#!/bin/bash

#check number of reads
zgrep -c "@" S436E_R1_trimmed.fastq.gz

#Save the headers from the S436 trimmed file and pipe it over to the gzip function to save as gzip file.
zgrep -h "@" S436E_R1_trimmed.fastq.gz | gzip -c > S436E_headers.fastq.gz

#double check to see if it ported it over correctly
zgrep -c "@" S436E_headers.fastq.gz > S436E_counts.txt
