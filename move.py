import os,re,sys
import gzip
import subprocess

## This script takes a list of file names, extracts the sample ID,
## and moves all matching files to a subdirectory for further analysis.

## format of files in list: unmapped_76_574_L01_FP100001024TR_1.fq
## sample ID: unmapped_76_574_L01_FP100001024TR

## Parent directory: /project/noujdine_61/mgosborn/Gametophytes/trimmed_reads_fastp

batchID = 'fourth'

batchName = batchID + 'Batch'
batchList = batchName + '.txt'

batchFolder = '/project/noujdine_61/mgosborn/Gametophytes/trimmed_reads_fastp/' + batchName

outNames = open('temp.txt','w')

with open(batchList, 'r') as files:
    for file in files:
        file_stripped = file.strip() #get rid of any trailing/leading white space
        file_name = file_stripped.rsplit('_',1)[0] #grab file name
        outNames.write(file_name + '\n')
    files.close()
outNames.close()

with open('temp.txt', 'r') as names:
    for name in names:
        name_stripped = name.strip()
        cmd = 'mv ' + name_stripped + '*.fq ' + batchFolder
        print(cmd)
        os.system(cmd)
    names.close()
