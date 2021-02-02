library(themetagenomics)
with(DAVID,table(META$Site,META$Donor))

CLEAN <- prepare_data(otu_table=DAVID$ABUND,
                      rows_are_taxa=FALSE,
                      tax_table=DAVID$TAX,
                      metadata=DAVID$META,
                      formula=~Site + s(Day),
                      refs='UBERON:saliva',
                      cn_normalize=FALSE,
                      drop=TRUE)

system.time(TOPICS <- find_topics(CLEAN,K=15))
system.time(TOPIC_EFFECTS <- est(TOPICS))
TOPIC_EFFECTS$topic_effects$`SiteUBERON:feces`$est
TOPIC_EFFECTS$topic_effects$Day$est
vis(TOPIC_EFFECTS,type='taxa')
vis(TOPIC_EFFECTS,type='continuous')

setwd("/Users/Mel/Desktop")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
library(microbiome)
library(ape)

otumat <- read.delim("AbundanceTable.txt", sep = "\t", row.names = 1, check.names = FALSE)
taxmat <- matrix(nrow = nrow(otumat), ncol = 0)
rownames(taxmat) <- rownames(otumat)
taxmat <- as.data.frame(cbind(Taxon = rownames(taxmat), taxmat))
levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
taxmat <- taxmat %>% separate(Taxon, levels, sep = ";")
taxmat <- as.matrix(taxmat)

#create phyloseq objects
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)

#Only keep Bacteria
physeq = subset_taxa(physeq, Domain=="Bacteria")

#### REMOVE SINGLETONS
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq) 

#### REMOVE UNCLASSIFIED PHYLUM
physeq <- subset_taxa(physeq, Phylum!="Unclassified Bacteria")

#load metadata csv file that denotes population
metadata = read.csv("10421_metadata.csv", row.names = 1)
metadata$Avg_Survival <- as.factor(as.character(metadata$Avg_Survival))

metadata_merged = select(metadata, -c("SampleName","Sequencer","Lane"))
metadata_merged = unique(metadata_merged)
rownames(metadata_merged) <- metadata_merged$SampleID

#add population info to physeq object
physeq@sam_data = sample_data(metadata)

#normalize abundance by sequencer
sequencers <- c(as.character(unique(metadata$Sequencer)))
for(sequencer in sequencers){
  temp = subset_samples(physeq, Sequencer==sequencer)
  
  total = median(sample_sums(temp))
  standf = function(x, t=total) round(t * (x / sum(x)))
  temp = transform_sample_counts(temp, standf)
  
  assign(sequencer,temp)
}

#list = cat(paste(shQuote(sequencers, type="cmd2"), collapse=", "))
physeq_merged = merge_phyloseq(FP100001023TR, FP100000947BR, FP100000945BR, FP100001024TR, FP100000946BR, FP100000944BR, FP100001026TR, FP100000948BR, FP100001025TR, V300048900, V300043045, V300049674, V300043035, DP8400010343BR, DP8400010332BR)

#Merge samples with multiple runs
physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
sample_data(physeq_merged) <- metadata_merged
#Only include planted individuals
physeq_merged = subset_samples(physeq_merged, Planted == "Planted") 

# center log transform
# physeq_merged_transformed <- microbiome::transform(physeq_merged, 'clr')
physeq_merged_transformed <- physeq_merged


### THEMETAGENOMICS

otumat <- as.data.frame(physeq_merged_transformed@otu_table)
taxmat <- as.data.frame(physeq_merged_transformed@tax_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
temp <- subset(metadata, select = "SampleID")
temp$rename <- rownames(temp)
rownames(otumat) <- temp$rename[match(rownames(otumat),temp$SampleID)]
rownames(metadata_merged) <- temp$rename[match(rownames(metadata_merged),temp$SampleID)]

## define a helper function
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

## transform all columns
taxmat <- taxmat %>% mutate_each(funs(empty_as_na)) 

CLEAN <- prepare_data(otu_table=otumat,
                      rows_are_taxa=FALSE,
                      tax_table=taxmat,
                      metadata=metadata_merged,
                      formula=~s(average_plant_biomass),
                      cn_normalize=FALSE,
                      drop=TRUE)

TOPICS <- find_topics(CLEAN,K=15)
TOPIC_EFFECTS <- est(TOPICS)

vis(TOPIC_EFFECTS, type = 'taxa')



