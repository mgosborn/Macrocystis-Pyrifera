import os,re,sys
import gzip
import subprocess

## This script runs kaiju for all files in a folder and merges the output in a consensus table.

## the files we want to run kaiju on are in the format:
## unmapped_sampleID_index_lane_readID_[1/2].fq
## examples:
### unmapped_9_510_L01_FP100001023TR_1.fq
### unmapped_9_510_L01_FP100001023TR_2.fq

print("Starting...")

print("Creating list of unmaped fqs in folder...")
cmd = "ls unmapped*_1.fq > files.txt"
os.system(cmd)

print("Grabbing file names...")

location_nodes_dmp = "/scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/nodes.dmp"
location_names_dmp = "/scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/names.dmp"
location_db_nr_euk_fmi = "/scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/kaiju_db_nr_euk.fmi"

with open('files.txt', 'r') as files:
    for file in files:
        file_stripped = file.strip() #get rid of any trailing/leading white space
        file_name = file_stripped.rsplit('_',1)[0] #grab file name
        print("Running kaiju on: " + file_name)
        cmd = "kaiju -z 25 -t " + location_nodes_dmp + ' -f ' + location_db_nr_euk_fmi +' -i ' + file_name + '_1.fq -j ' + file_name + '_2.fq -o kaiju_' + file_name + '.out'
        os.system(cmd)
        print("Making krona output...")
        cmd = "kaiju2krona -t " + location_nodes_dmp + " -n " + location_names_dmp + " -i kaiju_" + file_name + '.out -o kaiju_' + file_name + '.out.krona'
        os.system(cmd)
        print("Converting to html...")
        cmd = "ktImportText -o kaiju_" + file_name + '.out.html kaiju_' + file_name + '.out.krona'
        os.system(cmd)
    files.close()

#print("Individual kaiju complete. Merging profiles to a consensus table...")

#kaiju2table -t nodes.dmp -n names.dmp -p -o kaiju_summary.tsv kaiju.out [kaiju2.out, ...]



#kaiju -z 25 -t nodes.dmp -f kaiju_db_nr_euk.fmi -i inputfile.fastq -j secondfile.fastq -o kaiju.out
#kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.names.out
#kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona

#kaiju -z 25 -t /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/nodes.dmp -f /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/kaiju_db_nr_euk.fmi -i unmapped_123_529_L01_FP100000947BR_1.fq -j unmapped_123_529_L01_FP100000947BR_2.fq -o kaiju.out
#kaiju-addTaxonNames -t /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/nodes.dmp -n /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/names.dmp -i kaiju.out -o kaiju.names.out
#kaiju2krona -t /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/nodes.dmp -n /scratch/mgosborn/Gametophytes/PracticeSample/10samples/kaijudb/names.dmp -i kaiju.out -o kaiju.out.krona
