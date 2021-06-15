setwd("/Users/Mel/Desktop/Artemis/Kelp Microbiome")
library(phyloseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(microbiome)
library(vegan)
library(zetadiv)
library('vcfR')


## This script takes the abundance table and phenotype data and runs through entire microbiome/phenotype pipeline.

## Remember: Grab abundance table from 
## Command: metaxa2_dc -o AbundanceTable.txt -r "metaxa2_allreads" -p "^[^.]+" *level_7.txt

## Remember: Create updated phenotype data and plug in to line 40 below
## Use R script: CompileMetadataFromFarm.R


###################
###
### Loading abundance table and phenotype data. 
### Preparing phyloseq object.
###
###################

otumat <- read.delim("051721_AbundanceTable.txt", sep = "\t", row.names = 1, check.names = FALSE)
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

#### REMOVE SINGLETONS
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq) 

#load metadata csv file that denotes population
metadata = read.csv("051721_metadata.csv", row.names = 1)
metadata_merged = dplyr::select(metadata, -c("SampleName","Sequencer","Lane","Index"))
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
physeq_merged = merge_phyloseq(V300059128, V300062015, FP100000947BR, FP100001023TR, V300059139, V300059435, FP100000946BR, V300061963, V300059334, V300060570, V300062007, FP100000945BR, FP100000948BR, FP100001025TR, FP100001024TR, FP100001026TR, FP100000944BR, V300061032, V300059379, DP8400010332BR, V300059000, V300057490, V300058991, V300043045, V300043035, V300059290, V300058990, V300048900, DP8400010343BR, V300049674, V300058998, FP100001022TR)
rm(V300059128, V300062015, FP100000947BR, FP100001023TR, V300059139, V300059435, FP100000946BR, V300061963, V300059334, V300060570, V300062007, FP100000945BR, FP100000948BR, FP100001025TR, FP100001024TR, FP100001026TR, FP100000944BR, V300061032, V300059379, DP8400010332BR, V300059000, V300057490, V300058991, V300043045, V300043035, V300059290, V300058990, V300048900, DP8400010343BR, V300049674, V300058998, FP100001022TR)

#Merge samples with multiple runs
physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
# Since mean function is additive for abundance counts, divide by number of times sample occurs in dataset
temp_otu <- t(as.data.frame(physeq_merged@otu_table))
temp2 <- as.data.frame(table(metadata$SampleID))
temp2 <- subset(temp2, temp2$Freq>1)
row.names(temp2) <- temp2$Var1
for (z in temp2$Var1) {
  temp_otu[,z] <- round((temp_otu[,z])/temp2[z,"Freq"])
}
temp_otu <- t(temp_otu)

OTU = otu_table(temp_otu, taxa_are_rows = FALSE)
physeq_merged@otu_table = OTU
sample_data(physeq_merged) <- metadata_merged

#Only keep Bacteria
physeq_merged = subset_taxa(physeq_merged, Domain=="Bacteria")
physeq_merged = subset_taxa(physeq_merged, Phylum != "Unclassified Bacteria")
#topph = sort(tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Phylum"], sum), TRUE)[1:round(length(get_taxa_unique(physeq_merged, "Phylum"))/2)]
topph = sort(tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Phylum"], sum), TRUE)[1:5]
physeq_merged = subset_taxa(physeq_merged, Phylum %in% names(topph))


###################
###
### Barplots, NMDS, PCoA
###
###################

## BAR PLOT
#dev.new()
#plot_bar(physeq_merged, fill = "Phylum")

## PCA
temp <- microbiome::transform(physeq_merged, 'clr')
ord_clr <- phyloseq::ordinate(temp, "RDA", distance = "euclidean")

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
temp2 <- phyloseq::plot_ordination(physeq_merged, ord_clr, type="samples", color="Population", title = "PCA: Population (n = 559)") +
  geom_point(size = .5) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Population), linetype = 2)
dev.new()
print(temp2)

## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
# Generate distance matrix
temp_clr_dist_matrix <- phyloseq::distance(temp, method = "euclidean")

#ADONIS test
print("perMANOVA for Population (n = 559)")
print(vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Population))

# Check if ADONIS confounded by differences in dispersion
dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Population)
print(dispr)
print(permutest(dispr))

###################
###
### Subsetting kelp samples to top and bottom quartiles for subsequent analysis.
### Based on Total_Biomass_Rank.
###
###################

physeq_merged <- subset_samples(physeq_merged, Total_Biomass_Rank != "Middle 50%")

###################
###
### Kelp SNPs data (vcf format). Convert to genetic distance matrix for zeta diversity analysis.
###
###################

#load VCF
vcf <- read.vcfR("first_set_no_ghenghis_biallelic_pass_qc_max_missing_0.8_maf_0.02_meanmindp_5.vcf.gz")
#Convert VCF data to a genlight object
x <- vcfR2genlight(vcf)
#create distance matrix
x.dist <- dist(x, method = "euclidean")
x.dist.mat <- as.matrix(x.dist)

levels <- c("Phylum","Class","Order","Family","Genus","Species")

###################
###
### Pipeline: Alpha Div, Beta Div, Zeta Div, Permanova, Correlation, Diff. Abundance
###
###################

for(level in levels){
  #Use CLR to transform abundance table from counts to variance of each taxa relative to geometric mean of all taxa
  temp = tax_glom(physeq_merged, level)
  ##Use command below if want to investigate top half most abundant taxa
  #temp = prune_taxa(names(sort(taxa_sums(temp), TRUE))[1:round(ncol(temp@otu_table)/2)], temp)
  temp_no_trans <- temp # save for alpha div, zeta div, spieceasi, and aldex2 later on
  #subset for top and bottom quartiles based on biomass
  top25_no_trans <- subset_samples(temp_no_trans, Total_Biomass_Rank == sort(unique(temp_no_trans@sam_data[["Total_Biomass_Rank"]]))[2])
  bottom25_no_trans <- subset_samples(temp_no_trans, Total_Biomass_Rank == sort(unique(temp_no_trans@sam_data[["Total_Biomass_Rank"]]))[1])
  temp <- microbiome::transform(temp, 'clr')
  
  ###################
  ###
  ### ALPHA DIVERSITY
  ###
  ###################
  
  adiv <- data.frame(
    "Observed" = phyloseq::estimate_richness(temp_no_trans, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(temp_no_trans, measures = "Shannon"),
    "Chao1" = phyloseq::estimate_richness(temp_no_trans, measures = "Chao1"),
    "Total_Biomass_Rank" = phyloseq::sample_data(temp_no_trans)$Total_Biomass_Rank)

  #Plot adiv measures

  a <- adiv %>%
    gather(key = metric, value = value, c("Observed", "Shannon", "Chao1.Chao1")) %>%
    mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Chao1.Chao1"))) %>%
    ggplot(aes(x = Total_Biomass_Rank, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = Total_Biomass_Rank), height = 0, width = .2) +
    scale_color_manual(values=c("#FFC000", "#27CCBD")) +
    labs(x = "", y = "", title = paste("Taxa Level: ", level, " (Variable: Total_Biomass_Rank)",sep = "")) +
    facet_wrap(~ metric, scales = "free") +
    theme(legend.position="none")
  dev.new()
  print(a)
  
  ###################
  ###
  ### BETA DIVERSITY & PCA (Aitchison Distance)
  ###
  ###################
  
  # Aitchison Distance PCA Ordination (RDA without constraints is PCA)
  ord_clr <- phyloseq::ordinate(temp, "RDA", distance = "euclidean")
  ##Plot scree plot
  #dev.new()
  #phyloseq::plot_scree(ord_clr) +
  #  geom_bar(stat="identity", fill = "blue") +
  #  labs(x = "\nAxis", y = "Proportion of Variance\n", title = paste("Taxa Level: ",level,sep =""))
  
  ##Eigenvalues and % prop. variance explained
  #head(ord_clr$CA$eig)
  #sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))
  
  #Scale axes and plot ordination
  clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
  clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
  temp2 <- phyloseq::plot_ordination(temp, ord_clr, type="samples", color="Total_Biomass_Rank", title = paste("Taxa Level: ",level,sep ="")) +
    geom_point(size = 2) +
    scale_color_manual(values=c("#FFC000", "#27CCBD")) +
    coord_fixed(clr2 / clr1) +
    stat_ellipse(aes(group = Total_Biomass_Rank), linetype = 2)
  dev.new()
  print(temp2)
  
  ###################
  ###
  ### ZETA DIVERSITY 
  ###
  ###################
  
  ### Top 25%

  dev.new()
  z <- as.data.frame(top25_no_trans@otu_table)
  z[z >= 1] <- 1
  zeta_decline <- Zeta.decline.ex(data.spec = z, orders = 1:50, plot = TRUE)
  mtext(paste("Zeta Decline. Top 25% Biomass Quartile. Taxa Level: ", level, sep = ""), adj = 1.5, side = 3, line = 2.5)

  for (i in 2:5){
    dec <- Zeta.ddecay(xy = x.dist.mat, data.spec = z, order = i, distance.type = "custom", dist.custom = x.dist.mat, normalize = "Jaccard", plot = FALSE)
    dec <- dec$reg$model
    dec_name <- paste("dec",i,sep = "")
    assign(dec_name,dec)
  }

  dev.new()
  p <- ggplot() +
    stat_smooth(data = dec2, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec2")) +
    stat_smooth(data = dec3, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec3")) +
    stat_smooth(data = dec4, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec4")) +
    stat_smooth(data = dec5, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec5")) +
    labs(y = "Zeta diversity", x = "Genetic Distance (Euclidean)", title = paste("Zeta Decay. Top 25% Biomass Quartile. Taxa Level: ", level, sep = "")) +
    scale_color_discrete(name = "Zeta Order", labels = c("2", "3", "4","5"))
  print(p)

  ### Bottom 25%

  dev.new()
  z <- as.data.frame(bottom25_no_trans@otu_table)
  z[z >= 1] <- 1
  zeta_decline <- Zeta.decline.ex(data.spec = z, orders = 1:50, plot = TRUE)
  mtext(paste("Zeta Decline. Bottom 25% Biomass Quartile. Taxa Level: ", level, sep = ""), adj = 1.5, side = 3, line = 2.5)

  for (i in 2:5){
    dec <- Zeta.ddecay(xy = x.dist.mat, data.spec = z, order = i, distance.type = "custom", dist.custom = x.dist.mat, normalize = "Jaccard", plot = FALSE)
    dec <- dec$reg$model
    dec_name <- paste("dec",i,sep = "")
    assign(dec_name,dec)
  }

  dev.new()
  p <- ggplot() +
    stat_smooth(data = dec2, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec2")) +
    stat_smooth(data = dec3, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec3")) +
    stat_smooth(data = dec4, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec4")) +
    stat_smooth(data = dec5, method = "glm", aes(distance.reg, zeta.val.reg, color = "dec5")) +
    labs(y = "Zeta diversity", x = "Genetic Distance (Euclidean)", title = paste("Zeta Decay. Bottom 25% Biomass Quartile. Taxa Level: ", level, sep = "")) +
    scale_color_discrete(name = "Zeta Order", labels = c("2", "3", "4","5"))
  print(p)

  
  ###################
  ###
  ### PERMANOVA
  ###
  ###################
  
  ## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
  # Generate distance matrix
  temp_clr_dist_matrix <- phyloseq::distance(temp, method = "euclidean")

  #ADONIS test
  print(paste("perMANOVA for taxa level: ", level, sep = ""))
  print(vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Total_Biomass_Rank))

  # Check if ADONIS confounded by differences in dispersion
  dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Total_Biomass_Rank)
  print(dispr)
  #dev.new()
  #plot(dispr, main = paste("Ordination Centroids and Dispersion (Aitchison Distance). Taxa Level: ",level, sep = ""), sub = "")
  #boxplot(dispr, main = paste("Taxa Level: ", level, sep = ""), xlab = "")
  print(permutest(dispr))

  ###################
  ###
  ### CORRELATION (SpiecEasi)
  ###
  ###################
  
  #subset for top and bottom quartiles based on biomass
  top25 <- subset_samples(temp_no_trans, Total_Biomass_Rank == sort(unique(temp_no_trans@sam_data[["Total_Biomass_Rank"]]))[2])
  bottom25 <- subset_samples(temp_no_trans, Total_Biomass_Rank == sort(unique(temp_no_trans@sam_data[["Total_Biomass_Rank"]]))[1])
  
  #otumat <- as.data.frame(physeq_merged@otu_table)
  #namestop25 <- sample_names(top25)
  #namesbottom25 <- sample_names(bottom25)
  
  ####
  # top 25%
  ####
  
  se.mb.top25 <- spiec.easi(top25, method='mb')#,pulsar.params=list(thresh=0.01))
  
  optbeta <- as.matrix(symBeta(getOptBeta(se.mb.top25)))
  edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
  
  optbeta <- symBeta(getOptBeta(se.mb.top25))
  #edge_weights <- (1-Matrix::summary(t(optbeta))[,3])/2
  edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
  edge_weights <- edge_weights/max(edge_weights)
  
  #edge_weights <- (edge_weights/max(edge_weights))^2
  #edge_weights[edge_weights < summary(edge_weights)[5]] <- 0.05 ## zero all edge weights below third quartile
  
  
  ig2.mb.top25 <- adj2igraph(getRefit(se.mb.top25),  rmEmptyNodes=TRUE,
                             vertex.attr=list(name=taxa_names(top25)),
                             edge.attr=list(color=edge_cols, curved = .25, weight = edge_weights))
  
  df <- igraph::as_data_frame(ig2.mb.top25, 'both')
  check <- as.data.frame(taxmat[,"Phylum"])
  colnames(check) <- "Phylum"
  check$id <- row.names(check)
  
  df$vertices <- df$vertices %>% 
    left_join(check, c('name'='id'))
  
  updated_g <- graph_from_data_frame(df$edges,
                                     directed = F,
                                     vertices = df$vertices)
  #updated_g <- delete_edges(updated_g, which(E(updated_g)$weight == 0))
  #updated_g <- delete.vertices(simplify(updated_g), degree(updated_g)==0)
  
  
  updated_g$layout <- layout_with_dh
  deg <- degree(updated_g)
  
  coul <- scales::viridis_pal(option = "C")(length(unique(tax_table(top25)[,2])))
  
  edge_function <- paste(edge_cols,edge_weights,')',sep = "")
  edge_col_val <- c()
  for(i in 1:length(edge_function)){
    edge_col_val[i] <- eval(parse(text=edge_function[i]))
  }
  
  dev.new()
  #plot(updated_g, edge.width = E(updated_g)$weight, edge.color = edge_col_val, vertex.size = log10(deg)*6, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  #plot(updated_g, edge.color = edge_col_val, vertex.size = log10(deg)*6, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  plot(updated_g, edge.color = edge_col_val, vertex.size = log10(deg)*4, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  
  legend("bottomleft", legend=levels(as.factor(V(updated_g)$Phylum))  , col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 1 , horiz = FALSE, inset = c(-0.1, 0.2))
  
  
  ####
  # bottom 25%
  ####
  
  se.mb.bottom25 <- spiec.easi(bottom25, method='mb')#,pulsar.params=list(thresh=0.01))
  
  optbeta <- as.matrix(symBeta(getOptBeta(se.mb.bottom25)))
  edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]
  
  optbeta <- symBeta(getOptBeta(se.mb.bottom25))
  #edge_weights <- (1-Matrix::summary(t(optbeta))[,3])/2
  edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
  edge_weights <- edge_weights/max(edge_weights)
  
  #edge_weights <- (edge_weights/max(edge_weights))^2
  #edge_weights[edge_weights < summary(edge_weights)[5]] <- 0.05 ## zero all edge weights below third quartile
  
  
  ig2.mb.bottom25 <- adj2igraph(getRefit(se.mb.bottom25),  rmEmptyNodes=TRUE,
                                vertex.attr=list(name=taxa_names(bottom25)),
                                edge.attr=list(color=edge_cols, curved = .25, weight = edge_weights))
  
  df <- igraph::as_data_frame(ig2.mb.bottom25, 'both')
  check <- as.data.frame(taxmat[,"Phylum"])
  colnames(check) <- "Phylum"
  check$id <- row.names(check)
  
  df$vertices <- df$vertices %>% 
    left_join(check, c('name'='id'))
  
  updated_g <- graph_from_data_frame(df$edges,
                                     directed = F,
                                     vertices = df$vertices)
  #updated_g <- delete_edges(updated_g, which(E(updated_g)$weight == 0))
  #updated_g <- delete.vertices(simplify(updated_g), degree(updated_g)==0)
  
  
  updated_g$layout <- layout_with_dh
  deg <- degree(updated_g)
  
  coul <- scales::viridis_pal(option = "C")(length(unique(tax_table(bottom25)[,2])))
  
  edge_function <- paste(edge_cols,edge_weights,')',sep = "")
  edge_col_val <- c()
  for(i in 1:length(edge_function)){
    edge_col_val[i] <- eval(parse(text=edge_function[i]))
  }
  
  dev.new()
  #plot(updated_g, edge.width = E(updated_g)$weight, edge.color = edge_col_val, vertex.size = log10(deg)*6, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  #plot(updated_g, edge.color = edge_col_val, vertex.size = log10(deg)*6, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  plot(updated_g, edge.color = edge_col_val, vertex.size = log10(deg)*4, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))])
  
  legend("bottomleft", legend=levels(as.factor(V(updated_g)$Phylum))  , col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 1 , horiz = FALSE, inset = c(-0.1, 0.1))
  
  ###################
  ###
  ### DIFFERENTIAL ABUNDANCE (ALDEx2)
  ###
  ###################
  
  aldex2_da <- ALDEx2::aldex(t(data.frame(phyloseq::otu_table(temp_no_trans))), phyloseq::sample_data(temp_no_trans)$Total_Biomass_Rank, test="t", effect = TRUE, denom="iqlr")
  #dev.new()
  #ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

  #Clean up presentation
  sig_aldex2 <- aldex2_da %>%
    rownames_to_column(var = "OTU") %>%
    filter(wi.eBH < 0.05) %>%
    arrange(effect, wi.eBH) %>%
    dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
  print(sig_aldex2)
  
  ###################
  ###
  ### Save phyloseq object and rename
  ###
  ###################
  
  temp2 = paste("pm_clr_",level, sep = "")
  assign(temp2,temp)
  
}
