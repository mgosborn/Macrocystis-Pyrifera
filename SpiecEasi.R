setwd("/Volumes/Artemis/Kelp Microbiome/")
library("phyloseq"); packageVersion("phyloseq")
library(dplyr)
library(tidyr)
library(microbiome)
library(SpiecEasi)
library(igraph)

# ## Note: SpiecEasi pipeline transforms the data-- so input physeq object should NOT be clr transformed
# 
# ###################
# ###
# ### Loading abundance table and phenotype data. 
# ### Preparing phyloseq object.
# ###
# ###################
# 
# otumat <- read.delim("AbundanceTable.txt", sep = "\t", row.names = 1, check.names = FALSE)
# taxmat <- matrix(nrow = nrow(otumat), ncol = 0)
# rownames(taxmat) <- rownames(otumat)
# taxmat <- as.data.frame(cbind(Taxon = rownames(taxmat), taxmat))
# levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
# taxmat <- taxmat %>% separate(Taxon, levels, sep = ";")
# taxmat <- as.matrix(taxmat)
# 
# #create phyloseq objects
# OTU = otu_table(otumat, taxa_are_rows = TRUE)
# TAX = tax_table(taxmat)
# physeq = phyloseq(OTU, TAX)
# 
# #### REMOVE SINGLETONS
# physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq) 
# 
# #load metadata csv file that denotes population
# metadata = read.csv("20321_metadata.csv", row.names = 1)
# metadata$Avg_Survival <- as.factor(as.character(metadata$Avg_Survival))
# metadata_merged = dplyr::select(metadata, -c("SampleName","Sequencer","Lane"))
# metadata_merged = unique(metadata_merged)
# rownames(metadata_merged) <- metadata_merged$SampleID
# 
# #add population info to physeq object
# physeq@sam_data = sample_data(metadata)
# 
# #normalize abundance by sequencer
# sequencers <- c(as.character(unique(metadata$Sequencer)))
# for(sequencer in sequencers){
#   temp = subset_samples(physeq, Sequencer==sequencer)
#   
#   total = median(sample_sums(temp))
#   standf = function(x, t=total) round(t * (x / sum(x)))
#   temp = transform_sample_counts(temp, standf)
#   
#   assign(sequencer,temp)
# }
# 
# #list = cat(paste(shQuote(sequencers, type="cmd2"), collapse=", "))
# physeq_merged = merge_phyloseq(FP100001023TR, FP100000947BR, FP100001024TR, FP100000945BR, FP100000944BR, FP100000946BR, FP100001026TR, FP100001025TR, FP100000948BR, V300048900, V300043045, V300049674, DP8400010332BR, V300043035, DP8400010343BR, FP100001022TR)
# 
# #Merge samples with multiple runs
# physeq_merged = merge_samples(physeq_merged, "SampleID",fun = mean)
# # Since mean function is additive for abundance counts, divide by number of times sample occurs in dataset
# temp_otu <- t(as.data.frame(physeq_merged@otu_table))
# temp2 <- as.data.frame(table(metadata$SampleID))
# temp2 <- subset(temp2, temp2$Freq>1)
# row.names(temp2) <- temp2$Var1
# for (z in temp2$Var1) {
#   temp_otu[,z] <- round((temp_otu[,z])/temp2[z,"Freq"])
# }
# temp_otu <- t(temp_otu)
# 
# OTU = otu_table(temp_otu, taxa_are_rows = FALSE)
# physeq_merged@otu_table = OTU
# sample_data(physeq_merged) <- metadata_merged
# 
# #Only keep Bacteria
# physeq_merged = subset_taxa(physeq_merged, Domain=="Bacteria")
# 
# 



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
