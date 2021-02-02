library(dplyr)
library(tidyr)
library(tidyverse)

setwd("/Users/Mel/Desktop")

## This script takes the Code_SamplesDNA_for_SNPs.csv file (Gametophyte Code & Sample ID) 
## and Harvest_data.csv (raw phenotype data) and Harvest_data_survived_data_top_3_pheno_averaged.txt
## (processed phenotype data) and combines it with a list of sample names that analysis
## was done on. This becomes a metadata file that can be combined with a physeq object to make phyloseq graphs. 

## TO GENERATE LIST OF SAMPLE NAMES:
# Navigate to output directory /project/noujdine_61/mgosborn/Output
# List of files will be in format: metaxa2_unmapped_474_561_L01_FP100000945BR.level_7.txt
# type: ls *level_7* > names.txt

## Reformat list of file names to SampleID_Index_Lane_Sequencer
## i.e. metaxa2_unmapped_4_505_L01_FP100000946BR.level_7.txt to 4_505_L01_FP100000946BR
## Then split column to extract SampleID, which will be used to merge with other tables

names <- read.delim("names.txt", header = FALSE, stringsAsFactors = FALSE)
names <- as.data.frame(gsub("metaxa2_unmapped_","",names[,1]))
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
raw_pheno <- read.csv("Harvest_data.csv", header = TRUE)
proc_pheno <- read.delim("Harvest_data_survived_data_top_3_pheno_averaged.txt", header = TRUE, stringsAsFactors = FALSE)
proc_pheno <- proc_pheno[ , -which(names(proc_pheno) %in% c('X','X.1','blade','stipe','biomass','individual','Population'))]
names(proc_pheno)[names(proc_pheno)=="Gametophyte_code"] <- "GametophyteCode"

## Calculate phenotype averages and merge with raw_pheno dataframe

raw_pheno[ raw_pheno == "." ] <- NA
# 
# raw_pheno$Blade_weight <- as.numeric(as.character(raw_pheno$Blade_weight))
# raw_pheno$Stipe_weight <- as.numeric(as.character(raw_pheno$Stipe_weight))
# raw_pheno$Individual_count <- as.numeric(as.character(raw_pheno$Individual_count))
raw_pheno$Survival <- as.numeric(as.character(raw_pheno$Survival))
# 
# Avg_Blade_weight <- raw_pheno %>% 
#   group_by(GametophyteCode) %>% 
#   summarise(Avg_Blade_weight = mean(Blade_weight, na.rm = TRUE))
# 
# Avg_Stipe_Weight <- raw_pheno %>% 
#   group_by(GametophyteCode) %>% 
#   summarise(Avg_Stipe_Weight = mean(Stipe_weight, na.rm = TRUE))
# 
# Avg_Indiv_Count <- raw_pheno %>% 
#   group_by(GametophyteCode) %>% 
#   summarise(Avg_Indiv_Count = mean(Individual_count, na.rm = TRUE))
# 
Avg_Survival <- raw_pheno %>% 
  group_by(GametophyteCode) %>% 
  summarise(Avg_Survival = mean(Survival, na.rm = TRUE))

Avg_Survival <- as.data.frame(Avg_Survival)

#proc_pheno <- merge(proc_pheno, Avg_Survival, by )

# raw_pheno <- merge(raw_pheno, Avg_Blade_weight, by = "GametophyteCode")
# raw_pheno <- merge(raw_pheno, Avg_Stipe_Weight, by = "GametophyteCode")
# raw_pheno <- merge(raw_pheno, Avg_Indiv_Count, by = "GametophyteCode")
# raw_pheno <- merge(raw_pheno, Avg_Survival, by = "GametophyteCode")

## Reformat gametophyteID file
metadata <- merge(pops, names, by = "SampleID", no.dups = FALSE)
metadata$PopSex <- metadata$GametophyteCode
levels <- c("Population","temp","Sex","temp2")
metadata <- metadata %>% separate(PopSex, levels, sep = '\\.')
metadata <- metadata[ , -which(names(metadata) %in% c('temp','temp2'))]


#meta_inc_not_farmed <- metadata
#row.names(meta_inc_not_farmed) <- meta_inc_not_farmed$SampleName
#write.csv(meta_inc_not_farmed, "111320metadata.csv")
#names(proc_pheno)[names(proc_pheno)=="Gametophyte_code"] <- "GametophyteCode"
#test <- merge(proc_pheno, metadata, by = "GametophyteCode", all = TRUE)

#raw_pheno <- raw_pheno[ , -which(names(raw_pheno) %in% c("Population","Plant_ID","Replicate","Processing_date","Photo_ID", "Blade_weight","Stipe_weight","Individual_count","Survival","Sporophyte","Farm_code"))]
#raw_pheno <- raw_pheno[!duplicated(raw_pheno$GametophyteCode),]

metadata <- merge(metadata, proc_pheno, by = "GametophyteCode", all = TRUE)
metadata <- merge(metadata, Avg_Survival, by = "GametophyteCode", all = TRUE)
metadata <- metadata[complete.cases(metadata$SampleName),]

metadata$Planted <- ifelse(is.na(metadata$Avg_Survival), "Not Planted", "Planted")
row.names(metadata) <- metadata$SampleName

## Rank quantitative pheno values. Cutoff values are the average of the entire dataset (491 samples) calculated 11/16/20. 

metadata$Total_Biomass_Rank <- ifelse(metadata$average_plant_biomass > 168, "Big (> 168g)","Small (< 168g)")
metadata$Blade_Rank <- ifelse(metadata$average_blade_weight > 125, "Big Blades (> 125g)","Small Blades (< 125g)")
metadata$Stipe_Rank <- ifelse(metadata$average_Stipe_weight > 42, "Big Stipe (> 42g)", "Small Stipe (< 42g)")
metadata$Indiv_Rank <- ifelse(metadata$average_Individual_count > 4, "More individuals (> 4)","Less individuals (< 4)")

write.csv(metadata, "10721_metadata.csv")

