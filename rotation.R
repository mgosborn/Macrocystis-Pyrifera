setwd("/Users/Mel/Desktop")
library("phyloseq"); packageVersion("phyloseq")
library(psych)
library(microbiome)
library(vegan)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(Factoshiny)


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
metadata = read.csv("10721_metadata.csv", row.names = 1)
metadata$Avg_Survival <- as.factor(as.character(metadata$Avg_Survival))
metadata_merged = select(metadata, -c("SampleName","Sequencer","Lane"))
metadata_merged = unique(metadata_merged)
rownames(metadata_merged) <- metadata_merged$SampleID

#add population info to physeq object
physeq@sam_data = sample_data(metadata)
#####Only keep planted & survived individuals
physeq <- subset_samples(physeq, Planted == "Planted")
physeq <- subset_samples(physeq, is.na(average_plant_biomass) == FALSE) #keep who survived til harvest

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


##### PHYLOSEQ GRAPHS - pre transform
PCOA.ord <- ordinate(physeq_merged, "PCoA", "bray")
dev.new()
plot_ordination(physeq_merged, PCOA.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

NMDS.ord <- ordinate(physeq_merged, "NMDS", "bray") ##no convergence of stress ratio, not plotting
#dev.new()
#plot_ordination(physeq_merged, PCOA.ord, color = "Population", title = "Bray-Curtis NMDS: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

### CENTER LOG TRANSFORM 

physeq_merged_transformed <- microbiome::transform(physeq_merged, 'clr')
biomass <- as.data.frame(physeq_merged@sam_data)[,c("average_plant_biomass","average_Stipe_weight", "average_blade_weight")]
biomass_clr <- as.data.frame(microbiome::transform(biomass, 'clr'))


####### FACTOMINER

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(physeq_merged_transformed@otu_table)
# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,biomass_clr, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
#facto_mat <- na.omit(facto_mat)

### FactoMineR
res.pca <- PCA(facto_mat, quanti.sup=1272:1274, graph = FALSE)
res.shiny=PCAshiny(res.pca)
#invest <- Investigate(res.pca, document = "pdf_document")
#save <- dimdesc(res.pca, axes = 1:5)
######### Now do the same for only top 50% most abundant OTUs
PM_top50 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/2)], physeq_merged)
# center log transform
PM_top50_transformed <- microbiome::transform(PM_top50, 'clr')
# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top50_transformed@otu_table)
# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,biomass_clr, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
### FactoMineR
res.pca <- PCA(facto_mat, quanti.sup=637:639, graph = FALSE)
res.shiny=PCAshiny(res.pca)








# get otu, tax table, and metadata from physeq
spe_clr <- as.data.frame(physeq_merged_transformed@otu_table)

#parallel analysis to determine number of components to retain in PCA
pc_clr <- fa.parallel(spe_clr, fa="pc", n.iter=100,
                      show.legend=FALSE, main="CLR - Scree plot with parallel analysis")


######### Now do the same for only top 50% most abundant OTUs
PM_top50 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/2)], physeq_merged)
# center log transform
spe_clr_top50 <- microbiome::transform(PM_top50, 'clr')
# get otu, tax table, and metadata from physeq
spe_clr_top50 <- as.data.frame(spe_clr_top50@otu_table)
#parallel analysis to determine number of components to retain in PCA
pc_clr_top50 <- fa.parallel(spe_clr_top50, fa="pc", n.iter=100,
                      show.legend=FALSE, main="CLR Top 50%- Scree plot with parallel analysis")


######### Now do the same for only top 10% most abundant OTUs
PM_top10 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/10)], physeq_merged)
# center log transform
spe_clr_top10 <- microbiome::transform(PM_top10, 'clr')
# get otu, tax table, and metadata from physeq
spe_clr_top10 <- as.data.frame(spe_clr_top10@otu_table)
#parallel analysis to determine number of components to retain in PCA
pc_clr_top10 <- fa.parallel(spe_clr_top10, fa="pc", n.iter=100,
                            show.legend=FALSE, main="CLR Top 10%- Scree plot with parallel analysis")

######### Now do the same for only top 5% most abundant OTUs
PM_top5 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/20)], physeq_merged)
# center log transform
spe_clr_top5 <- microbiome::transform(PM_top5, 'clr')
# get otu, tax table, and metadata from physeq
spe_clr_top5 <- as.data.frame(spe_clr_top5@otu_table)
#parallel analysis to determine number of components to retain in PCA
pc_clr_top5 <- fa.parallel(spe_clr_top5, fa="pc", n.iter=100,
                            show.legend=FALSE, main="CLR Top 5%- Scree plot with parallel analysis")


######### Now do the same for only top 1% most abundant OTUs
PM_top1 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/100)], physeq_merged)
# center log transform
spe_clr_top1 <- microbiome::transform(PM_top1, 'clr')
# get otu, tax table, and metadata from physeq
spe_clr_top1 <- as.data.frame(spe_clr_top1@otu_table)
#parallel analysis to determine number of components to retain in PCA
pc_clr_top1 <- fa.parallel(spe_clr_top1, fa="pc", n.iter=100,
                            show.legend=FALSE, main="CLR Top 1%- Scree plot with parallel analysis")



### Extracting Principal components
pc_all_none <- principal(spe_clr, nfactors=36, rotate="none", scores=TRUE)
pc_all_varimax <- principal(spe_clr, nfactors=36, rotate="varimax", scores=TRUE)
pc_all_promax <- principal(spe_clr, nfactors=36, rotate="promax", scores=TRUE)
temp <- cor(biomass_clr$average_plant_biomass, pc_all_none$score)

dev.new()
factor.plot(pc_all_none)

pc_50_none <- principal(spe_clr_top50, nfactors=28, rotate="none", scores=TRUE)
pc_50_varimax <- principal(spe_clr_top50, nfactors=28, rotate="varimax", scores=TRUE)
pc_50_promax <- principal(spe_clr_top50, nfactors=28, rotate="promax", scores=TRUE)
cor(biomass_clr$average_plant_biomass, pc_50_promax$score)








# combine otu and metadata - prepare for FactoMineR
#facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
#row.names(facto_mat) <- facto_mat$Row.names
#facto_mat <- subset(facto_mat, select = -c(Row.names))
#facto_mat <- na.omit(facto_mat)
#explore <- subset(facto_mat, select = -c(Population))

## CLR applied to biomass data? 
#facto_mat[1353:1355] <- microbiome::transform(facto_mat[1353:1355], 'clr')