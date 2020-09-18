import os,re,sys
import gzip
import subprocess

## This script runs metaphlan for all files in a folder and merges the output in a consensus table.
## This script must be run in the metaPhlan-env.

## the files we want to run metaphlan on are in the format:
## unmapped_sampleID_index_lane_readID_[1/2].fq
## examples:
### unmapped_9_510_L01_FP100001023TR_1.fq
### unmapped_9_510_L01_FP100001023TR_2.fq

print("Starting...")

print("Creating list of unmaped fqs in folder...")
cmd = "ls unmapped* > files.txt"
print(cmd)
os.system(cmd)

print("Grabbing file names...")

with open('files.txt', 'r') as files:
    for file in files:
        file_stripped = file.strip() #get rid of any trailing/leading white space
        file_name = file_stripped.split('.fq')[0] #grab file name
        print("Running metaphlan on new file...")
        cmd = "metaphlan " + file_stripped + ' --bowtie2out ' + file_name +'.bowtie2.bz2 --nproc 5 --input_type fastq --unknown_estimation -o profiled_' + file_name + '.txt'
        os.system(cmd)
        
    files.close()

print("Individual metaphlan complete. Merging profiles to a single abundance table...")

open('merged_abundance_table.txt','w')

cmd = "merge_metaphlan_tables.py profiled*.txt > merged_abundance_table.txt"
os.system(cmd)
