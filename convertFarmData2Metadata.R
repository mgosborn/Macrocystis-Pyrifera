library(dplyr)
library(tidyr)

setwd("/Users/Mel/Desktop")

## This script takes the Code_SamplesDNA_for_SNPs.csv file and combines it with
## a list of sample names that analysis was done on. This becomes a metadata file
## that is consequently loaded as a phyloseq object to make phyloseq graphs. 

## TO GENERATE LIST OF SAMPLE NAMES:
# Navigate to output directory /project/noujdine_61/mgosborn/Output
# List of files will be in format: metaxa2_unmapped_474_561_L01_FP100000945BR.level_8.txt
# type: ls *level_8* > names.txt

## GametophyteCode and SampleID file
## GametophyteCode format: POP.#.F.alpha#
## i.e. AQ.02.F.B2
## SampleID matches SampleID in names file
pops <- read.csv("Code_SamplesDNA_for_SNPs.csv", header = TRUE)

##List of file names in the format SampleNumber_Index_Lane_Sequencer
## i.e. 4_505_L01_FP100000946BR
names <- read.delim("names.txt", header = FALSE, stringsAsFactors = FALSE)
names <- as.data.frame(gsub("metaxa2_unmapped_","",names[,1]))
names <- as.data.frame(gsub(".level_8.txt","",names[,1]))

## Reformat names file
names$Names <- names[,1] #duplicate column
levels <- c("SampleID","Index","Lane","Sequencer")
names <- names %>% separate(Names, levels, sep = "_") #Split Names by '_'
colnames(names)[1] <- "SampleName" #rename first column

## Reformat gametophyteID file
pops <- read.csv("Code_SamplesDNA_for_SNPs.csv", header = TRUE)
pops$Population <- pops[,1]
pops$Population <- gsub("\\...*","",pops$Population) #get rid of everything after first '.'

metadata <- merge(pops, names, by = "SampleID")
metadata$SampleID <- as.numeric(as.character(metadata$SampleID))
metadata <- metadata[order(metadata$SampleID),] #Sort matrix by SampleID column (ascending)
row.names(metadata) <- metadata$SampleName
write.csv(metadata, "110520metadata.csv")
