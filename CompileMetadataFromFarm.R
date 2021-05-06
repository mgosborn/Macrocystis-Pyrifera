library(dplyr)
library(tidyr)
library(tidyverse)

setwd("/Users/Mel/Desktop")

## This script takes four input files: 
## (1) Code_SamplesDNA_for_SNPs.csv (Gametophyte Code & Sample ID), 
## (2) Harvest_data.csv (raw biomass data)
## (3) Plant_Chl.csv (raw chlorophyll data)
## (4) names.txt (described below)
## ... to generate a single metadata file that can be combined with a phyloseq object for future analysis.

## TO GENERATE LIST OF SAMPLE NAMES FROM METAXA2 OUTPUT (names.txt):
# Navigate to output directory
# First Set: /project/noujdine_61/mgosborn/Gametophytes/trimmed_reads_fastp/renamed_allreads_fa/level_7_output
# Second Set: TBD
# List of files for the first set will be in format: metaxa2_allreads_101_507_L01_FP100000947BR.level_7.txt
# type: ls *level_7* > allreads_names.txt

## Reformat list of file names to SampleID_Index_Lane_Sequencer
## i.e. metaxa2_allreads_101_507_L01_FP100000947BR to 101_507_L01_FP100000947BR
## Then split column to extract SampleID, which will be used to merge with other tables

names <- read.delim("allreads_names.txt", header = FALSE, stringsAsFactors = FALSE)
names <- as.data.frame(gsub("metaxa2_allreads_","",names[,1]))
names <- as.data.frame(gsub(".level_7.txt","",names[,1]))
names$Names <- names[,1] #duplicate column
levels <- c("SampleID","Index","Lane","Sequencer")
names <- names %>% separate(Names, levels, sep = "_") #Split Names by '_'
colnames(names)[1] <- "SampleName" #rename first column

## GametophyteCode and SampleID file
## GametophyteCode format: POP.#.[FM].alpha#
## i.e. AQ.02.F.B2
## SampleID matches SampleID in names file

pops <- read.csv("Code_SamplesDNA_for_SNPs.csv", header = TRUE)

pheno <- read.csv("Harvest_data.csv", header = TRUE, stringsAsFactors = FALSE)
#names(pheno)[names(pheno)=="Farm_code"] <- "SampleID"
pheno[ pheno == "." ] <- NA
temp <- c("Blade_weight","Stipe_weight")
pheno[temp] <- sapply(pheno[temp],as.numeric)
pheno$Total_Biomass <- pheno$Blade_weight + pheno$Stipe_weight

avg <- pheno %>% # average the top 3 values (if available, otherwise the average of what is available)
  group_by(Farm_code) %>% 
  top_n(n = 3, wt = Blade_weight) %>%
  summarize(average_blade_weight=mean(Blade_weight))
pheno <- merge(pheno,avg, by = "Farm_code", all = TRUE)

avg <- pheno %>% # average the top 3 values (if available, otherwise the average of what is available)
  group_by(Farm_code) %>% 
  top_n(n = 3, wt = Stipe_weight) %>%
  summarize(average_stipe_weight=mean(Stipe_weight))
pheno <- merge(pheno,avg, by = "Farm_code", all = TRUE)

avg <- pheno %>% # average the top 3 values (if available, otherwise the average of what is available)
  group_by(Farm_code) %>% 
  top_n(n = 3, wt = Total_Biomass) %>%
  summarize(average_total_biomass=mean(Total_Biomass))
pheno <- merge(pheno,avg, by = "Farm_code", all = TRUE)
rm(avg)

#rm(pheno)

## Chlorophyll A data
chla <- read.csv("Plant_Chl.csv", header = TRUE)
names(chla)[names(chla)=="Plant"] <- "Plant_ID"
chla <- chla[!grepl("WHITE|BLACK|BAD", chla$Plant_ID),] #remove bad or missing data
chla <- chla[!grepl("NaN", chla$ChlA),] #remove bad or missing data

levels <- c("Farm_code","Replicate")
chla <- chla %>% separate(Plant_ID, levels, sep = '(?<=[0-9])(?=[A-Z])') #separate sample ID (#) and replicate (alpha)
chla <- na.omit(chla) #remove missing data
chla <- chla[!as.numeric(chla$Farm_code) > 500,] # remove samples with known typos... 

chla <- chla %>% # average the top 3 values (if available, otherwise the average of what is available)
  group_by(Farm_code) %>% 
  top_n(n = 3) %>%
  summarize(average_chla=mean(ChlA))

pheno <- merge(pheno, chla, by = "Farm_code", all = TRUE)
proc_pheno <- pheno[,-which(names(pheno) %in% c('Farm_code','Plant_ID','Replicate','Processing_date','Blade_weight','Stipe_weight','Photo_ID','Individual_count','Survival','Total_Biomass'))]
proc_pheno <- distinct(proc_pheno)

## Compile metadata 
metadata <- merge(pops,proc_pheno, by = "GametophyteCode", all = TRUE)
metadata <- merge(names,metadata, by = "SampleID", all = TRUE)
metadata <- metadata[complete.cases(metadata$SampleName),]
metadata$Population <- ifelse(is.na(metadata$Population), sub("\\..*", "", metadata$GametophyteCode), metadata$Population)
metadata$Planted <- ifelse(is.na(metadata$average_blade_weight), "Not Planted", "Planted")
row.names(metadata) <- metadata$SampleName


## Rank quantitative phenotype values.  

metadata$Total_Biomass_Rank <- ifelse(metadata$average_total_biomass > summary(metadata$average_total_biomass)[5], paste("Top 25% (> ",round(summary(metadata$average_total_biomass)[5],2),"g)", sep = ""),
                                      ifelse(metadata$average_total_biomass < summary(metadata$average_total_biomass)[2], paste("Bottom 25% (< ",round(summary(metadata$average_total_biomass)[2],2),"g)", sep = ""), "Middle 50%"))
metadata$Blade_Rank <- ifelse(metadata$average_blade_weight > summary(metadata$average_blade_weight)[5], paste("Top 25% (> ",round(summary(metadata$average_blade_weight)[5],2),"g)", sep = ""),
                              ifelse(metadata$average_blade_weight < summary(metadata$average_blade_weight)[2], paste("Bottom 25% (< ",round(summary(metadata$average_blade_weight)[2],2),"g)", sep = ""), "Middle 50%"))
metadata$Stipe_Rank <- ifelse(metadata$average_stipe_weight > summary(metadata$average_stipe_weight)[5], paste("Top 25% (> ",round(summary(metadata$average_stipe_weight)[5],2),"g)", sep = ""),
                              ifelse(metadata$average_stipe_weight < summary(metadata$average_stipe_weight)[2], paste("Bottom 25% (< ",round(summary(metadata$average_stipe_weight)[2],2),"g)", sep = ""), "Middle 50%"))
metadata$Chla_Rank <- ifelse(metadata$average_chla > summary(metadata$average_chla)[5], paste("Top 25% (> ",round(summary(metadata$average_chla)[5],2),")", sep = ""),
                              ifelse(metadata$average_chla < summary(metadata$average_chla)[2], paste("Bottom 25% (< ",round(summary(metadata$average_chla)[2],2),")", sep = ""), "Middle 50%"))


write.csv(metadata, "050621_metadata.csv")
