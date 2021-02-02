setwd("/Users/Mel/Desktop")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(dplyr)
library(tidyr)
library(microbiome)
library(ape)


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

# center log transform
physeq_merged_transformed <- microbiome::transform(physeq_merged, 'clr')
#physeq_merged_transformed <- transform_sample_counts(physeq_merged_transformed, physeq_merged_transformed[physeq_merged_transformed < 0.0] <- 0.0)
physeq_merged_transformed@otu_table[physeq_merged_transformed@otu_table < 0.0] <- 0.0

pmt_genus = tax_glom(physeq_merged_transformed, "Genus")
pmt_family = tax_glom(physeq_merged_transformed, "Family")
pmt_phylum = tax_glom(physeq_merged_transformed, "Phylum")

PCOA.ord <- ordinate(pmt_phylum, "PCoA", "bray")
dev.new()
plot_ordination(pmt_phylum, PCOA.ord, color = "Population", title = "Bray-Curtis NMDS: Population - Phylum") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 


pmt_phylum_planted = subset_samples(pmt_phylum, Planted == "Planted")
pmt_genus_planted = subset_samples(pmt_genus, Planted == "Planted")
pmt_family_planted = subset_samples(pmt_family, Planted == "Planted")

PCOA.ord <- ordinate(pmt_family_planted, "PCoA", "bray")
dev.new()
plot_ordination(pmt_family_planted, PCOA.ord, color = "average_plant_biomass", title = "Bray-Curtis PCoA: Avg Biomass - Family") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")


## Keep only the most abundant phyla (top 50%).
phyla.sum = tapply(taxa_sums(pmt_phylum_planted), tax_table(pmt_phylum_planted)[, "Phylum"], sum, na.rm=TRUE)
topHalfPhyla = names(sort(phyla.sum, TRUE))[1:round(ncol(pmt_phylum_planted@otu_table)/10)]
PMtopHalf = prune_taxa((tax_table(pmt_phylum_planted)[, "Phylum"] %in% topHalfPhyla), pmt_phylum_planted)

PCoAPMtopHalf.ord <- ordinate(PMtopHalf, "PCoA", "bray")
dev.new()
plot_ordination(PMtopHalf, PCoAPMtopHalf.ord, color = "Population", title = "Bray-Curtis PCoA: Phyla Top Tenth") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 
dev.new()
plot_ordination(PMtopHalf, PCoAPMtopHalf.ord, color = "average_plant_biomass", title = "Bray-Curtis PCoA: Phyla Top Tenth") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")





#normalized bar plot without black bars between unique taxa
dev.new()
plot_bar(physeq_merged_transformed, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

PCoA.ord <- ordinate(physeq_merged_transformed, "PCoA", "bray")
plot_ordination(physeq_merged_transformed, PCoA.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

physeq_mt_PLANTED_ONLY = subset_samples(physeq_merged_transformed, Planted == "Planted")

PCoA_planted.ord <- ordinate(physeq_mt_PLANTED_ONLY, "PCoA", "bray")
dev.new()
plot_ordination(physeq_mt_PLANTED_ONLY, PCoA_planted.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

dev.new()
plot_ordination(physeq_mt_PLANTED_ONLY, PCoA_planted.ord, color = "average_plant_biomass", title = "Bray-Curtis PCoA: Avg Biomass") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")


## Keep only the most abundant species (top 50%).
species.sum = tapply(taxa_sums(physeq_mt_PLANTED_ONLY), tax_table(physeq_mt_PLANTED_ONLY)[, "Species"], sum, na.rm=TRUE)
topHalfSpecies = names(sort(species.sum, TRUE))[1:round(ncol(physeq_mt_PLANTED_ONLY@otu_table)/2)]
PMtopHalf = prune_taxa((tax_table(physeq_mt_PLANTED_ONLY)[, "Species"] %in% topHalfSpecies), physeq_mt_PLANTED_ONLY)

PCoAPMtopHalf.ord <- ordinate(PMtopHalf, "PCoA", "bray")
plot_ordination(PMtopHalf, PCoAPMtopHalf.ord, color = "Population", title = "Bray-Curtis PCoA: Population Top Half") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 


## Keep only the most abundant species (top 25%).
species.sum = tapply(taxa_sums(physeq_mt_PLANTED_ONLY), tax_table(physeq_mt_PLANTED_ONLY)[, "Species"], sum, na.rm=TRUE)
topQuarterSpecies = names(sort(species.sum, TRUE))[1:round(ncol(physeq_mt_PLANTED_ONLY@otu_table)/4)]
PMtop25 = prune_taxa((tax_table(physeq_mt_PLANTED_ONLY)[, "Species"] %in% topHalfSpecies), physeq_mt_PLANTED_ONLY)

PCoAPMtop25.ord <- ordinate(PMtop25, "PCoA", "bray")
dev.new()
plot_ordination(PMtop25, PCoAPMtop25.ord, color = "Population", title = "Bray-Curtis PCoA: Population Top Quarter") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 


## Keep only the most abundant species (top 10%).

species.sum = tapply(taxa_sums(physeq_mt_PLANTED_ONLY), tax_table(physeq_mt_PLANTED_ONLY)[, "Species"], sum, na.rm=TRUE)
topTenthSpecies = names(sort(species.sum, TRUE))[1:round(ncol(physeq_mt_PLANTED_ONLY@otu_table)/10)]
PMtopTenth = prune_taxa((tax_table(physeq_mt_PLANTED_ONLY)[, "Species"] %in% topTenthSpecies), physeq_mt_PLANTED_ONLY)

PCoAPMtopTenth.ord <- ordinate(PMtopTenth, "PCoA", "bray")
dev.new()
plot_ordination(PMtopTenth, PCoAPMtopTenth.ord, color = "Population", title = "Bray-Curtis PCoA: Population Top Tenth") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

dev.new()
plot_ordination(PMtopTenth, PCoAPMtopTenth.ord, color = "average_plant_biomass", title = "Bray-Curtis PCoA: Avg Biomass Top Tenth") + geom_point(size=2)

dev.new()
plot_ordination(PMtopTenth, PCoAPMtopTenth.ord, color = "average_blade_weight", title = "Bray-Curtis PCoA: Avg Blade Top Tenth") + geom_point(size=2)

dev.new()
plot_ordination(PMtopTenth, PCoAPMtopTenth.ord, color = "average_Stipe_weight", title = "Bray-Curtis PCoA: Avg Stipe Top Tenth") + geom_point(size=2)

## Keep only the most abundant species (top 5%).

species.sum = tapply(taxa_sums(physeq_mt_PLANTED_ONLY), tax_table(physeq_mt_PLANTED_ONLY)[, "Species"], sum, na.rm=TRUE)
topFiveSpecies = names(sort(species.sum, TRUE))[1:round(ncol(physeq_mt_PLANTED_ONLY@otu_table)/20)]
PMtopFive = prune_taxa((tax_table(physeq_mt_PLANTED_ONLY)[, "Species"] %in% topFiveSpecies), physeq_mt_PLANTED_ONLY)
PCoAPMtopFive.ord <- ordinate(PMtopFive, "PCoA", "bray")
dev.new()
plot_ordination(PMtopFive, PCoAPMtopFive.ord, color = "Population", title = "Bray-Curtis PCoA: Population Top Five") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

dev.new()
plot_ordination(PMtopFive, PCoAPMtopFive.ord, color = "average_plant_biomass", title = "Bray-Curtis PCoA: Avg Biomass Top Five") + geom_point(size=2)

dev.new()
plot_ordination(PMtopFive, PCoAPMtopFive.ord, color = "average_blade_weight", title = "Bray-Curtis PCoA: Avg Blade Top Five") + geom_point(size=2)

dev.new()
plot_ordination(PMtopFive, PCoAPMtopFive.ord, color = "average_Stipe_weight", title = "Bray-Curtis PCoA: Avg Stipe Top Five") + geom_point(size=2)






dev.new()
plot_heatmap(physeq_merged_transformed, taxa.label="Phylum")

PCoAphyseq.ord <- ordinate(physeq_merged, "PCoA", "bray")

plot_ordination(physeq_merged, PCoAphyseq.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 




physeq_mt_PLANTED_ONLY = subset_samples(physeq_merged_transformed, Planted == "Planted")

PCoAphyseq_mt.ord <- ordinate(physeq_merged_transformed, "PCoA", "bray")
PCoAphyseq_mt_PLANTED.ord <- ordinate(physeq_PLANTED_ONLY, "PCoA", "bray")

plot_ordination(physeq_merged_transformed, PCoAphyseq.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED.ord, color = "average_plant_biomass", shape = "Population", title = "Bray-Curtis PCoA: Population & Total Biomass") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED.ord, color = "Total_Biomass_Rank", shape = "Population", title = "Bray-Curtis PCoA: Population & Total Biomass") + geom_point(size=2)






#normalize abundance
total = median(sample_sums(pruned_physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
pruned_physeq = transform_sample_counts(pruned_physeq, standf)

NMDSphyseq.ord <- ordinate(pruned_physeq, "NMDS", "bray")
PCoAphyseq.ord <- ordinate(pruned_physeq, "PCoA", "bray")

#NMDS/PCoA Ordination - Population
dev.new()
plot_ordination(physeq, PCoAphyseq.ord, color = "Population", title = "Bray-Curtis PCoA: Population")




physeq_PLANTED_ONLY = subset_samples(physeq, Planted == "Planted")




#Calculate bray-curtis NMDS and PCoA ordination 
NMDSphyseq.ord <- ordinate(physeq, "NMDS", "bray")
NMDSphyseq_PLANTED_ONLY.ord <- ordinate(physeq_PLANTED_ONLY, "NMDS", "bray")
PCoAphyseq.ord <- ordinate(physeq, "PCoA", "bray")
PCoAphyseq_PLANTED_ONLY.ord <- ordinate(physeq_PLANTED_ONLY, "PCoA", "bray")

#NMDS/PCoA Ordination - Population
dev.new()
plot_ordination(physeq, NMDSphyseq.ord, color = "Population", title = "Bray-Curtis NMDS: Population")
dev.new()
plot_ordination(physeq, PCoAphyseq.ord, color = "Population", title = "Bray-Curtis PCoA: Population") + geom_point(size=2) + scale_color_manual(values = c("darkslategray3","darksalmon","gray40","mediumorchid1")) 

#NMDS/PCoA Ordination - Population and Sex
dev.new()
plot_ordination(physeq, NMDSphyseq.ord, color = "Sex", shape = "Population", title = "Bray-Curtis NMDS: Population & Sex")
plot_ordination(physeq, PCoAphyseq.ord, color = "Sex", shape = "Population", title = "Bray-Curtis PCoA: Population & Sex") + geom_point(size=2)

#NMDS Ordination - Population and Total Biomass
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "average_plant_biomass", shape = "Population", title = "Bray-Curtis NMDS: Population & Total Biomass")
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "Total_Biomass_Rank", shape = "Population", title = "Bray-Curtis NMDS: Population & Total Biomass")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "average_plant_biomass", shape = "Population", title = "Bray-Curtis PCoA: Population & Total Biomass") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "Total_Biomass_Rank", shape = "Population", title = "Bray-Curtis PCoA: Population & Total Biomass") + geom_point(size=2)

#NMDS Ordination - Population and Blade Weight
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "average_blade_weight", shape = "Population", title = "Bray-Curtis NMDS: Population & Blade Weight")
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "Blade_Rank", shape = "Population", title = "Bray-Curtis NMDS: Population & Blade Weight")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "average_blade_weight", shape = "Population", title = "Bray-Curtis PCoA: Population & Blade Weight") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "Blade_Rank", shape = "Population", title = "Bray-Curtis PCoA: Population & Blade Weight") + geom_point(size=2)

#NMDS Ordination - Population and Stipe Weight
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "average_Stipe_weight", shape = "Population", title = "Bray-Curtis NMDS: Population & Stipe Weight")
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "Stipe_Rank", shape = "Population", title = "Bray-Curtis NMDS: Population & Stipe Weight")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "average_Stipe_weight", shape = "Population", title = "Bray-Curtis PCoA: Population & Stipe Weight") + geom_point(size=2) + scale_color_gradient(low = "blue", high = "orange")
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "Stipe_Rank", shape = "Population", title = "Bray-Curtis PCoA: Population & Stipe Weight") + geom_point(size=2)

#NMDS Ordination - Population and Survival Rate
dev.new()
plot_ordination(physeq_PLANTED_ONLY, NMDSphyseq_PLANTED_ONLY.ord, color = "Avg_Survival", shape = "Population", title = "Bray-Curtis NMDS: Population & Survival Rate") + scale_color_manual(values = c("slategray","slategray3", "turquoise4","turquoise3","violetred4","violetred1")) 
plot_ordination(physeq_PLANTED_ONLY, PCoAphyseq_PLANTED_ONLY.ord, color = "Avg_Survival", shape = "Population", title = "Bray-Curtis PCoA: Population & Survival Rate") + geom_point(size=2) + scale_color_manual(values = c("slategray","slategray3", "turquoise4","turquoise3","violetred4","violetred1")) 








#### REMOVE SINGLETONS

pruned_physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 
pruned_physeq_PLANTED_ONLY <- prune_taxa(taxa_sums(physeq_PLANTED_ONLY) > 2, physeq_PLANTED_ONLY) 


#Calculate bray-curtis NMDS and PCoA ordination 
pruned_NMDSphyseq.ord <- ordinate(pruned_physeq, "NMDS", "bray")
pruned_NMDSphyseq_PLANTED_ONLY.ord <- ordinate(pruned_physeq_PLANTED_ONLY, "NMDS", "bray")
pruned_PCoAphyseq.ord <- ordinate(pruned_physeq, "PCoA", "bray")
pruned_PCoAphyseq_PLANTED_ONLY.ord <- ordinate(pruned_physeq_PLANTED_ONLY, "PCoA", "bray")

#NMDS/PCoA Ordination - Population
dev.new()
plot_ordination(physeq, NMDSphyseq.ord, color = "Population", title = "Bray-Curtis NMDS: Population (No 10-tons)")
dev.new()
plot_ordination(physeq, PCoAphyseq.ord, color = "Population", title = "Bray-Curtis PCoA: Population (No 10-tons)")









plot_ordination(physeq, physeq.ord, type = "samples", color = "Population") + geom_point(size=2) + scale_color_manual(values = c("chartreuse4", "thistle4","dodgerblue3")) 

plot_heatmap(physeq, taxa.label="Phylum", sample.label = "Population")






#normalized heat map
dev.new()
plot_heatmap(physeq, taxa.label="Phylum")

#species richness
dev.new()
plot_richness(physeq, measures = c("Observed","Shannon"))

