import os,re,sys
import gzip
import subprocess

print("Starting...")

## This script renames and converts a list of fq files into fa files, then moves fa files into separate folders for subsequent analysis
## fq files of interest here are the original (unfiltered) M. pyrifera gametophyte WGS reads

#Input Location: /project/noujdine_61/mgosborn/Gametophytes/trimmed_reads_fastp
#Output Location: '/project/noujdine_61/mgosborn/Gametophytes/trimmed_reads_fastp/renamed_allreads_fa'

## Format of list_of_sample_and_sequence_file_info.txt

# ./sample_id/read_lane_indexNum_fwdFileSuffix
# ./9/FP100000947BR_L01_510_1.fq.gz

## Format of fq files

# read_lane_indexNum_fwdFileSuffix
# read_lane_indexNum_revFileSuffix
# FP100000947BR_L01_510_1.fq.gz
# FP100000947BR_L01_510_2.fq.gz

ListOfReeads = open('list_of_sample_and_sequence_file_info.txt') #list of files that you want to rename and move
fwdFileSuffix = '1.fq.gz'
revFileSuffix = '2.fq.gz'

for line in ListOfReeads: #for each fq filename in list of reads
    print(line)
    line_stripped = line.strip() # get rid of any trailing/leading white space
    sample_id = line_stripped.split('/')[1] #grab sample id
    #print(sample_id)
    raw_fastq_filename = line_stripped.split('/')[2]
    prefix_of_fastq = raw_fastq_filename.split('1.fq.gz')[0]
    #print(prefix_of_fastq)
    FWDfq = str(prefix_of_fastq) + str(fwdFileSuffix)
    REVfq = str(prefix_of_fastq) + str(revFileSuffix)
    #print(FWDfq)
    #print(REVfq)
    read = raw_fastq_filename.split('_')[0]
    lane = raw_fastq_filename.split('_')[1]
    indexNum = raw_fastq_filename.split('_')[2]
    #print(read, lane, indexNum)
    #index_lane_read = str(indexNum) + '_' + str(lane) + '_' + str(read)

    FWDallreadsNewFilename = 'allreads_' + str(sample_id) + '_' + str(indexNum) + '_' + str(lane) + '_' + str(read) + '_1.fa'
    REVallreadsNewFilename = 'allreads_' + str(sample_id) + '_' + str(indexNum) + '_' + str(lane) + '_' + str(read) + '_2.fa'

    with open('list_of_renamed_allreads_fa.txt', 'a') as files: #keep track of new files
        files.write(FWDallreadsNewFilename)
        files.close()
    
    print("Converting FWD fq to fa (seqtk)...")
    cmd = "seqtk seq -a " + FWDfq + " > " + FWDallreadsNewFilename
    os.system(cmd)

    print("Converting REV fq to fa (seqtk)...")
    cmd = "seqtk seq -a " + REVfq + " > " + REVallreadsNewFilename
    os.system(cmd)

    print("Done. Next:")

print("Finished renaming and converting to fa. Moving to new folder called renamed_allreads_fa")

cmd = "mkdir renamed_allreads_fa"
os.system(cmd)

cmd = "mv allreads*fa > renamed_allreads_fa"
os.system(cmd)

cmd = "mv list_of_renamed_allreads_fa.txt > renamed_allreads_fa"
os.system(cmd)
