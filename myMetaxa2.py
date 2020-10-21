import os,re,sys
import gzip
import subprocess

## This script runs metaxa2 for all files in a folder and merges the output in a consensus table.

## the files we want to run metaxa2 on are in the format:
## unmapped_sampleID_index_lane_readID_[1/2].fq
## examples:
### unmapped_9_510_L01_FP100001023TR_1.fq
### unmapped_9_510_L01_FP100001023TR_2.fq

#metaxa2 -1 [first] -2 [second] -o [outputname] -f p --cpu 12 --plus T

TaxaLevel = "6"
# Kingdom/Domain (1), Phylum (2), Class (3), Order (4), Family (5), Genus (6) and Species (7)

print("Starting...")

print("Creating list of unmaped fqs in folder...")
cmd = "ls unmapped*_1.fq > files.txt"
os.system(cmd)

print("Grabbing file names...")

with open('files.txt', 'r') as files:
    for file in files:
        file_stripped = file.strip() #get rid of any trailing/leading white space
        file_name = file_stripped.rsplit('_',1)[0] #grab file name
        metaxa2output_name = "metaxa2_" + file_name
        print("Running metaxa2 on: " + file_name)
        cmd = "metaxa2 -1 " + file_name + "_1.fq -2 " + file_name + "_2.fq -o " + metaxa2output_name + " -f p --cpu 12 --plus T"
        print(cmd)
        os.system(cmd)
        print("Counting taxa in " + file_name)
        cmd = "metaxa2_ttt -i " + metaxa2output_name + ".taxonomy.txt -o " + metaxa2output_name
        os.system(cmd)
    files.close()

print("Creating consensus abundance table...")
cmd = "metaxa2_dc -o AbundanceTable.txt -p \"^...\" *level_" + TaxaLevel + ".txt"
print(cmd)
os.system(cmd)


