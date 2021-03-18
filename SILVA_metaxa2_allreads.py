import os,re,sys
import gzip
import subprocess

## This script runs metaxa2 for all files in a folder and merges the output in a consensus table.

## the files we want to run metaxa2 on are in the format:
## allreads_sampleID_index_lane_readID_[1/2].fa
## examples:
### allreads_9_510_L01_FP100001023TR_1.fa
### allreads_9_510_L01_FP100001023TR_2.fa

#metaxa2 -1 [first] -2 [second] -o [outputname] -f pa --cpu 12 --plus T

# TaxaLevel = "7"
# Kingdom/Domain (1), Phylum (2), Class (3), Order (4), Family (5), Genus (6) and Species (7)

print("Starting...")

print("Removing temporary files leftover from halted metaxa2 run...")
cmd = "rm -rf metaxa_temp_*"
os.system(cmd)

print("Creating list of allreads fas in folder...")
cmd = "ls allreads*_1.fa > files.txt"
os.system(cmd)

print("Creating list of samples already done...")
cmd = "ls *level_7* > temp.txt"
os.system(cmd)

with open('temp.txt','r') as temps:
    with open('done.txt','w') as dones:
        for temp in temps:
            temp_stripped = temp.strip()
            temp_name = temp_stripped.rsplit('.level',1)[0]
            temp_temp = temp_name.strip("metaxa2_")
            temp_matchname = temp_temp + '_1.fa\n'
            dones.write(temp_matchname)
        dones.close()
    temps.close()
os.remove("temp.txt")

print("Ignoring files already done...")

cmd = "grep -xvf done.txt files.txt > remaining.txt"
os.system(cmd)

print("Made list of remaining files to run metaxa2 on.")

with open('remaining.txt', 'r') as files:
    for file in files:
        file_stripped = file.strip() #get rid of any trailing/leading white space
        file_name = file_stripped.rsplit('_',1)[0] #grab file name
        metaxa2output_name = "metaxa2_" + file_name
        print("Running metaxa2 on: " + file_name)
        cmd = "metaxa2 -g SSU_SILVA128 -1 " + file_name + "_1.fa -2 " + file_name + "_2.fa -o " + metaxa2output_name + " -f a --cpu 12 --plus T"
        print(cmd)
        os.system(cmd)
        print("Counting taxa in " + file_name)
        cmd = "metaxa2_ttt -i " + metaxa2output_name + ".taxonomy.txt -o " + metaxa2output_name
        print(cmd)
        os.system(cmd)
    files.close()

#print("Creating consensus abundance table...")
## COMMAND: metaxa2_dc -o AbundanceTable.txt -r "metaxa2_allreads_" -p "^[^.]+" *level_7.txt
## string after "-p" indicates sample name in matrix should be everything before the first dot in the filename" 
#cmd = "metaxa2_dc -o AbundanceTable.txt -r \"metaxa2_allreads_\" -p \"^[^.]+\" *level_" + TaxaLevel + ".txt"
#print(cmd)
#os.system(cmd)
