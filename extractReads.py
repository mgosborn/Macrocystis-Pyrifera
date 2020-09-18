import os,re,sys
import gzip
import subprocess

print("Starting...")

## This script extracts a lsit of unmapped reads from a list of corresponding fq files.

## Formats:

# ./sample_id/read_lane_indexNum_fwdFileSuffix
# ./9/FP100000947BR_L01_510_1.fq.gz

# read_lane_indexNum_fwdFileSuffix
# read_lane_indexNum_revFileSuffix
# FP100000947BR_L01_510_1.fq.gz
# FP100000947BR_L01_510_2.fq.gz

# lane_index_read_UnmappedTrail
# 510_L01_FP100000947BR_CI_03_hi_c_200404.unmapped.readnames

ListOfReeads = open('10samples.txt')
UnmappedTrail = '_CI_03_hi_c_200404.unmapped.readnames'
fwdFileSuffix = '1.fq.gz'
revFileSuffix = '2.fq.gz'
newFQunmapped = 'unmapped_'
UnmappedReadsFolder = '/scratch/mgosborn/Gametophytes/unmapped_readnames_for_melisa/'

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
    index_lane_read = str(indexNum) + '_' + str(lane) + '_' + str(read)
    unmapped_readname_file = str(index_lane_read) + UnmappedTrail
    unmapped_read_location = str(UnmappedReadsFolder) + unmapped_readname_file
    #print(unmapped_readname_file)
    # Make new name for new file that will contain unmapped reads only
    FWDunmappedNewFilename = 'unmapped_' + str(sample_id) + '_' + str(indexNum) + '_' + str(lane) + '_' + str(read) + '_1.fq'
    REVunmappedNewFilename = 'unmapped_' + str(sample_id) + '_' + str(indexNum) + '_' + str(lane) + '_' + str(read) + '_2.fq'

    print("Making new temp files...")
    FWDnewread = ""
    REVnewread = ""
    with open(unmapped_read_location, 'r') as reads:
        for read in reads:
            FWDnewread+=read.strip()+"/1\n"
            REVnewread+=read.strip()+"/2\n"
        reads.close()
    with open('FWD_temp_unmapped.txt', 'w') as reads:
        reads.write(FWDnewread)
        reads.close()
    with open('REV_temp_unmapped.txt','w') as reads:
        reads.write(REVnewread)
        reads.close()

    FWD_temp_unmapped = 'FWD_temp_unmapped.txt'
    REV_temp_unmapped = 'REV_temp_unmapped.txt'
    
    print("Doing FWD seqtk...")
    cmd = "seqtk subseq " + FWDfq + " " + FWD_temp_unmapped + " > " + FWDunmappedNewFilename
    os.system(cmd)

    print("Doing REV seqtk...")
    cmd = "seqtk subseq " + REVfq + " " + REV_temp_unmapped + " > " + REVunmappedNewFilename
    os.system(cmd)
    
    print("Done. Next:")
   
