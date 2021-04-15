setwd("/Volumes/Artemis/Kelp Microbiome/03_03_21 Gametophyte Analysis")
library('vcfR')
library('adegenet')
library("phyloseq"); packageVersion("phyloseq")
library(dplyr)
library(tidyr)
library(microbiome)
library(zetadiv)


####################
#### 
#### Abundance table
#### 
####################

### PRE-PROCESS 

## REMEMBER: Get most up to date Abundance Table. 
## COMMAND: metaxa2_dc -o AbundanceTable.txt -r "metaxa2_unmapped_" -p "^[^.]+" *level_7.txt

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

#### REMOVE SINGLETONS & DOUBLETONS
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 

#load metadata csv file that denotes population
metadata = read.csv("20321_metadata.csv", row.names = 1)
metadata$Avg_Survival <- as.factor(as.character(metadata$Avg_Survival))
metadata_merged = dplyr::select(metadata, -c("SampleName","Sequencer","Lane"))
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
physeq_merged = merge_phyloseq(FP100001023TR, FP100000947BR, FP100001024TR, FP100000945BR, FP100000944BR, FP100000946BR, FP100001026TR, FP100001025TR, FP100000948BR, V300048900, V300043045, V300049674, DP8400010332BR, V300043035, DP8400010343BR, FP100001022TR)

#Merge samples with multiple runs
physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
sample_data(physeq_merged) <- metadata_merged

#### REMOVE SINGLETONS & DOUBLETONS
physeq_merged <- prune_taxa(taxa_sums(physeq_merged) > 2, physeq_merged) 
top25 <- subset_samples(physeq_merged, Total_Biomass_Rank == "Top 25% (> 215.17g)")
bottom25 <- subset_samples(physeq_merged, Total_Biomass_Rank == "Bottom 25% (< 62.17g)")

temp <- tax_glom(top25, taxrank = "Species", NArm = TRUE)
species_top25 <- as.data.frame(temp@otu_table)
temp <- tax_glom(bottom25, taxrank = "Species", NArm = TRUE)
species_bottom25 <- as.data.frame(temp@otu_table)

temp <- tax_glom(top25, taxrank = "Genus", NArm = TRUE)
genus_top25 <- as.data.frame(temp@otu_table)
temp <- tax_glom(bottom25, taxrank = "Genus", NArm = TRUE)
genus_bottom25 <- as.data.frame(temp@otu_table)

temp <- tax_glom(top25, taxrank = "Phylum", NArm = TRUE)
phylum_top25 <- as.data.frame(temp@otu_table)
temp <- tax_glom(bottom25, taxrank = "Phylum", NArm = TRUE)
phylum_bottom25 <- as.data.frame(temp@otu_table)

#otumat <- as.data.frame(physeq_merged@otu_table)
namestop25 <- sample_names(top25)
namesbottom25 <- sample_names(bottom25)



####################
#### 
#### Creating the distance matrix
#### 
####################

#load VCF 
vcf <- read.vcfR("first_set_no_ghenghis_biallelic_pass_qc_max_missing_0.8_maf_0.02_meanmindp_5.vcf.gz")
top25vcf <- read.vcfR("TOP25only.vcf.gz")
bottom25vcf <- read.vcfR("BOTTOM25only.vcf.gz")

#Convert VCF data to a genlight object
x <- vcfR2genlight(vcf)
x.genind <- vcfR2genind(vcf)
x.genclone <- poppr::as.genclone(x.genind)

#create distance matrix
x.dist <- dist(x, method = "euclidean")
x.dist.mat <- as.matrix(x.dist)
x.nei.dist <- poppr::nei.dist(x.genclone, warning = TRUE)

x.dist <- diss.dist(x.genind)

####################
#### 
#### Zeta diversity decay
#### 
####################
speciesmat
genusmat
phylummat
x.dist.mat

#convert to presence/absence counts i.e. all numbers >= 1 become 1
species_top25[species_top25 >= 1] <- 1
species_bottom25[species_bottom25 >= 1] <- 1
genus_top25[genus_top25 >= 1] <- 1
genus_bottom25[genus_bottom25 >= 1] <- 1
phylum_top25[phylum_top25 >= 1] <- 1
phylum_bottom25[phylum_bottom25 >= 1] <- 1


dev.new()
species_decline <- Zeta.decline.ex(data.spec = speciesmat, orders = 1:50, plot = TRUE)
genus_decline <- Zeta.decline.ex(data.spec = genusmat, orders = 1:50, plot = TRUE)
phylum_decline <- Zeta.decline.ex(data.spec = phylummat, orders = 1:50, plot = TRUE)

dev.new()
species_t25_decline <- Zeta.decline.ex(data.spec = species_top25, orders = 1:50, plot = TRUE)
genus_t25_decline <- Zeta.decline.ex(data.spec = genus_top25, orders = 1:50, plot = TRUE)
phylum_t25_decline <- Zeta.decline.ex(data.spec = phylum_top25, orders = 1:50, plot = TRUE)

dev.new()
species_b25_decline <- Zeta.decline.ex(data.spec = species_bottom25, orders = 1:50, plot = TRUE)
genus_b25_decline <- Zeta.decline.ex(data.spec = genus_bottom25, orders = 1:50, plot = TRUE)
phylum_b25_decline <- Zeta.decline.ex(data.spec = phylum_bottom25, orders = 1:50, plot = TRUE)

dev.new()
species_decay <- Zeta.ddecays(xy = x.dist.mat, data.spec = speciesmat, orders = 2:15, sam = 50000, distance.type = "custom", dist.custom = x.dist.mat)
genus_decay <- Zeta.ddecays(xy = x.dist.mat, data.spec = genusmat, orders = 2:15, sam = 50000, distance.type = "custom", dist.custom = x.dist.mat)
phylum_decay <- Zeta.ddecays(xy = x.dist.mat, data.spec = phylummat, orders = 2:15, sam = 50000, distance.type = "custom", dist.custom = x.dist.mat)

dev.new()
species_dec2 <- Zeta.ddecay(xy = x.dist.mat, data.spec = speciesmat, order = 2, distance.type = "custom", dist.custom = x.dist.mat)
species_dec3 <- Zeta.ddecay(xy = x.dist.mat, data.spec = speciesmat, order = 3, distance.type = "custom", dist.custom = x.dist.mat)
species_dec4 <- Zeta.ddecay(xy = x.dist.mat, data.spec = speciesmat, order = 4, distance.type = "custom", dist.custom = x.dist.mat)
species_dec5 <- Zeta.ddecay(xy = x.dist.mat, data.spec = speciesmat, order = 5, distance.type = "custom", dist.custom = x.dist.mat)

dev.new()
genus_dec2 <- Zeta.ddecay(xy = x.dist.mat, data.spec = genusmat, order = 2, distance.type = "custom", dist.custom = x.dist.mat)
genus_dec3 <- Zeta.ddecay(xy = x.dist.mat, data.spec = genusmat, order = 3, distance.type = "custom", dist.custom = x.dist.mat)
genus_dec4 <- Zeta.ddecay(xy = x.dist.mat, data.spec = genusmat, order = 4, distance.type = "custom", dist.custom = x.dist.mat)
genus_dec5 <- Zeta.ddecay(xy = x.dist.mat, data.spec = genusmat, order = 5, distance.type = "custom", dist.custom = x.dist.mat)

dev.new()
phylum_dec2 <- Zeta.ddecay(xy = x.dist.mat, data.spec = phylummat, order = 2, distance.type = "custom", dist.custom = x.dist.mat)
phylum_dec3 <- Zeta.ddecay(xy = x.dist.mat, data.spec = phylummat, order = 3, distance.type = "custom", dist.custom = x.dist.mat)
phylum_dec4 <- Zeta.ddecay(xy = x.dist.mat, data.spec = phylummat, order = 4, distance.type = "custom", dist.custom = x.dist.mat)
phylum_dec5 <- Zeta.ddecay(xy = x.dist.mat, data.spec = phylummat, order = 5, distance.type = "custom", dist.custom = x.dist.mat)
