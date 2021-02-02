setwd("/Users/Mel/Desktop")
library("phyloseq"); packageVersion("phyloseq")
library(microbiome)
library(vegan)
library(FactoMineR)
library(Factoshiny)

#library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
#library(ape)



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

#### REMOVE UNCLASSIFIED PHYLUM
#physeq <- subset_taxa(physeq, Phylum!="Unclassified Bacteria")

#load metadata csv file that denotes population
metadata = read.csv("10721_metadata.csv", row.names = 1)
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
physeq_merged = merge_phyloseq(FP100001023TR, FP100000947BR, FP100001024TR, FP100000945BR, FP100000944BR, FP100000946BR, FP100001026TR, FP100001025TR, FP100000948BR, V300048900, V300043045, V300049674, DP8400010332BR, V300043035, DP8400010343BR, FP100001022TR)

#Merge samples with multiple runs
physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
sample_data(physeq_merged) <- metadata_merged

#### REMOVE SINGLETONS & DOUBLETONS
physeq_merged <- prune_taxa(taxa_sums(physeq_merged) > 2, physeq_merged) 

### CENTER LOG TRANSFORM 

# center log transform
physeq_merged_transformed <- microbiome::transform(physeq_merged, 'clr')
#physeq_merged_transformed <- transform_sample_counts(physeq_merged_transformed, physeq_merged_transformed[physeq_merged_transformed < 0.0] <- 0.0)
#physeq_merged_transformed@otu_table[physeq_merged_transformed@otu_table < 0.0] <- 0.0

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(physeq_merged_transformed@otu_table)
#taxmat <- as.data.frame(physeq_merged_transformed@tax_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)
#explore <- subset(facto_mat, select = -c(Population))

## CLR applied to biomass data? 
#facto_mat[1353:1355] <- microbiome::transform(facto_mat[1353:1355], 'clr')


### FactoMineR
res.pca <- PCA(facto_mat, quanti.sup=1353:1355, quali.sup=1352, graph = FALSE)
res.shiny=PCAshiny(res.pca)
#invest <- Investigate(res.pca, document = "pdf_document")
save <- dimdesc(res.pca, axes = 1:5)


######### Now do the same for only top 50% most abundant OTUs
PM_top50 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/2)], physeq_merged)

# center log transform
PM_top50_transformed <- microbiome::transform(PM_top50, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top50_transformed@otu_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

## CLR applied to biomass data? 
#facto_mat[678:680] <- microbiome::transform(facto_mat[678:680], 'clr')

### FactoMineR
PM_top50.pca <- PCA(facto_mat, quanti.sup=678:680, quali.sup=677, graph = FALSE)
PM_top50.shiny=PCAshiny(PM_top50.pca)
save_PM_top50 <- dimdesc(PM_top50.pca, axes = 1:5, proba = 0.1)




#########  Top 10% of OTUs
PM_top10 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/10)], physeq_merged)
# center log transform
PM_top10_transformed <- microbiome::transform(PM_top10, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top10_transformed@otu_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

## CLR applied to biomass data? 
#facto_mat[137:139] <- microbiome::transform(facto_mat[137:139], 'clr')

### FactoMineR
PM_top10.pca <- PCA(facto_mat, quanti.sup=137:139, quali.sup=136, graph = FALSE)
PM_top10.shiny=PCAshiny(PM_top10.pca)
save_PM_top10 <- dimdesc(PM_top10.pca, axes = 1:5, proba = 0.1)

#########  Top 5% of OTUs
PM_top5 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/20)], physeq_merged)

# center log transform
PM_top5_transformed <- microbiome::transform(PM_top5, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top5_transformed@otu_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

## CLR applied to biomass data? 
facto_mat[70:72] <- microbiome::transform(facto_mat[70:72], 'clr')

### FactoMineR
PM_top5.pca <- PCA(facto_mat, quanti.sup=70:72, quali.sup=69, graph = FALSE)
PM_top5.shiny=PCAshiny(PM_top5.pca)
save_PM_top5 <- dimdesc(PM_top5.pca, axes = 1:5, proba = 0.1)

#########  Top 1% of OTUs
PM_top1 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/100)], physeq_merged)

# center log transform
PM_top1_transformed <- microbiome::transform(PM_top1, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top1_transformed@otu_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

## CLR applied to biomass data? 
facto_mat[16:18] <- microbiome::transform(facto_mat[16:18], 'clr')

### FactoMineR
PM_top1.pca <- PCA(facto_mat, quanti.sup=16:18, quali.sup=15, graph = FALSE)
PM_top1.shiny=PCAshiny(PM_top1.pca)
save_PM_top1 <- dimdesc(PM_top1.pca, axes = 1:5, proba = 0.1)

#figure out top 10 contributors to variance from 

head1 <- as.matrix(head(save$Dim.1$quanti, n = 10), rownames = TRUE)
tail1 <- as.matrix(tail(save$Dim.1$quanti, n = 10), rownames = TRUE)
top1 <- abs(as.data.frame(rbind(head1,tail1)))
top1 <- top1[order(-top1$correlation),]
top1 <- top1[1:10,]

head2 <- as.matrix(head(save$Dim.2$quanti, n = 10), rownames = TRUE)
tail2 <- as.matrix(tail(save$Dim.2$quanti, n = 10), rownames = TRUE)
top2 <- abs(as.data.frame(rbind(head2,tail2)))
top2 <- top2[order(-top2$correlation),]
top2 <- top2[1:10,]

head3 <- as.matrix(head(save$Dim.3$quanti, n = 10), rownames = TRUE)
tail3 <- as.matrix(tail(save$Dim.3$quanti, n = 10), rownames = TRUE)
top3 <- abs(as.data.frame(rbind(head3,tail3)))
top3 <- top3[order(-top3$correlation),]
top3 <- top3[1:10,]

head4 <- as.matrix(head(save$Dim.4$quanti, n = 10), rownames = TRUE)
tail4 <- as.matrix(tail(save$Dim.4$quanti, n = 10), rownames = TRUE)
top4 <- abs(as.data.frame(rbind(head4,tail4)))
top4 <- top4[order(-top4$correlation),]
top4 <- top4[1:10,]

head5 <- as.matrix(head(save$Dim.5$quanti, n = 10), rownames = TRUE)
tail5 <- as.matrix(tail(save$Dim.5$quanti, n = 10), rownames = TRUE)
top5 <- abs(as.data.frame(rbind(head5,tail5)))
top5 <- top5[order(-top5$correlation),]
top5 <- top5[1:10,]

top_taxa <- row.names(top1)
top_taxa <- append(top_taxa,row.names(top2))
top_taxa <- append(top_taxa,row.names(top3))
top_taxa <- append(top_taxa,row.names(top4))
top_taxa <- append(top_taxa,row.names(top5))



pruned <- prune_taxa(top_taxa,physeq_merged)
temp <- prune_taxa(taxa_sums(pruned) > 1, pruned)
temp = filter_taxa(pruned,mean(taxa_sums(pruned)) > 1e-5, TRUE)
#physeq_merged <- prune_taxa(taxa_sums(physeq_merged) > 4, physeq_merged) 

# center log transform
pruned_transformed <- microbiome::transform(pruned, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(pruned_transformed@otu_table)
#taxmat <- as.data.frame(pruned_transformed@tax_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)
#explore <- subset(facto_mat, select = -c(Population))

### FactoMineR
res.pca <- PCA(facto_mat, quanti.sup=32:34, quali.sup=31, graph = FALSE)
res.shiny=PCAshiny(res.pca)
save2 <- dimdesc(res.pca, axes = 1:5)

#interactive pca
#res.shiny=PCAshiny(res.pca)

dev.new()
plot.PCA(res.pca,choix='var',col.quanti.sup='#0000FF')
#plot.PCA(res.pca,invisible=c('ind.sup'),label =c('ind','quali'))























#### TRY WITH TOP 5% ABUNDANT GENERA


## Keep only the most abundant genera (top 5%).
phyla.sum = tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Genus"], sum, na.rm=TRUE)
topTenPhyla = names(sort(phyla.sum, TRUE))[1:round(ncol(physeq_merged@otu_table)/20)]
PMtopTen = prune_taxa((tax_table(physeq_merged)[, "Genus"] %in% topTenPhyla), physeq_merged)

PMtopTen_transformed <- microbiome::transform(PMtopTen, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PMtopTen_transformed@otu_table)
taxmat <- as.data.frame(PMtopTen_transformed@tax_table)
metadata_merged<-metadata_merged[metadata_merged$Planted=="Planted",]
sub_metadata <- subset(metadata_merged, select = c("Population","average_plant_biomass","average_Stipe_weight","average_blade_weight"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)
explore <- subset(facto_mat, select = -c(Population))

### FactoMineR
res.pca <- PCA(facto_mat, quanti.sup=574:576, quali.sup=573)
res.shiny=PCAshiny(res.pca)









### AITCHISON DISTANCE 

# Aitchison is the euclidean of clr data 

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
veg_otu <- psotu2veg(physeq_merged_transformed)

# Calculate distance matrix
veg_otu_distmat <- vegdist(veg_otu, method = "euclidean")
aitch_distmat <- as.matrix(veg_otu_distmat)


