library(phyloseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(microbiome)
library(vegan)
library(zetadiv)
library(SpiecEasi)
library(igraph)



## This script runs through network and zeta analysis. 

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

#load metadata csv file that denotes population
metadata = read.csv("070721_metadata.csv", row.names = 1)
metadata_merged = dplyr::select(metadata, -c("SampleName","Sequencer","Lane","Index"))
metadata_merged = unique(metadata_merged)
rownames(metadata_merged) <- metadata_merged$SampleID

#add population info to physeq object
physeq@sam_data = sample_data(metadata)

#### Remove singletons and doubletons
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 

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

#Merge samples with multiple runs (warnings will be produced, don't worry- it is regarding the metadata. will be fixed in next lines)
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

#### Remove singletons and doubletons
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq)

#update OTU table and metadata in physeq object
OTU = otu_table(temp_otu, taxa_are_rows = FALSE)
physeq_merged@otu_table = OTU
sample_data(physeq_merged) <- metadata_merged

#Only keep  bacteria
physeq_merged = subset_taxa(physeq_merged, Domain=="Bacteria")
#physeq_merged = subset_taxa(physeq_merged, Phylum != "Unclassified Bacteria")

#Only keep top 10 most abundant Phyla. 
topph = sort(tapply(taxa_sums(physeq_merged), tax_table(physeq_merged)[, "Phylum"], sum), TRUE)[1:10]
physeq_merged = subset_taxa(physeq_merged, Phylum %in% names(topph))

###################
###
### Creating separate phyloseq objects for all populations and all quartiles. 
### Total of 20 phyloseq objects. 
###
###################

pops = c("All", "LC", "CB", "CI", "AQ")
levels <- c("Class","Order","Family","Genus","Species")
zeta_order = 10
z_decline = data.frame()
z_decay = data.frame()
x.dist.mat = read.csv("euclidean_dist_combined_SNPs.csv")
fraction_unclassified = data.frame()
pos_and_neg_edges = data.frame()
hubscores = data.frame()

for (p in pops){
  if (p == "All"){
    physeq_subpop = physeq_merged
  } else
    physeq_subpop = subset_samples(physeq_merged, Population == p)
  
  for (level in levels){
    #Conglomerate to specified taxa level
    physeq_subpop_glom = tax_glom(physeq_subpop, level)
    #Get rid of taxa unclassified at that level
    unclassified = str_subset(as.data.frame(physeq_subpop_glom@tax_table)[[level]], "Unclassified.*")
    physeq_subpop_glom = prune_taxa(!tax_table(physeq_subpop_glom)[,level] %in% unclassified, physeq_subpop_glom)
    
    temp = c(p, level, (length(unclassified)/(nrow(physeq_subpop_glom@tax_table)+length(unclassified))))
    fraction_unclassified <- rbind(temp, fraction_unclassified)
    
    # Create separate physeq objects for each quartile (based on total biomass)
    quantiles = quantile(physeq_subpop_glom@sam_data[["average_total_biomass"]], na.rm = TRUE)
    physeq_subpop_glom_1 <- subset_samples(physeq_subpop_glom, average_total_biomass <= quantiles[2])
    physeq_subpop_glom_2 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[2] & average_total_biomass <= quantiles[3])
    physeq_subpop_glom_3 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[3] & average_total_biomass <= quantiles[4])
    physeq_subpop_glom_4 <- subset_samples(physeq_subpop_glom, average_total_biomass > quantiles[4])
    
    #Add metadata
    physeq_subpop_glom_1@sam_data$Quantile = paste("Quantile 1 ( <", round(quantiles[2],2), "g)")
    physeq_subpop_glom_2@sam_data$Quantile = paste("Quantile 2 (", round(quantiles[2],2), "g to ",round(quantiles[3], 2), "g)")
    physeq_subpop_glom_3@sam_data$Quantile = paste("Quantile 3 (", round(quantiles[3],2), "g to ",round(quantiles[4], 2), "g)")
    physeq_subpop_glom_4@sam_data$Quantile = paste("Quantile 4 ( >", round(quantiles[4],2), "g)")
    
    quantiles = c("physeq_subpop_glom_1","physeq_subpop_glom_2","physeq_subpop_glom_3","physeq_subpop_glom_4")
    
    #Combine into single physeq object
    physeq_subpop_glom <- merge_phyloseq(physeq_subpop_glom_1, physeq_subpop_glom_2, physeq_subpop_glom_3, physeq_subpop_glom_4)
    
    ###################
    ###
    ### ALPHA DIVERSITY
    ###
    ###################
    
    if (p == pops[1]){
      adiv <- data.frame(
        "Observed Richness" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Observed"),
        "Population" = phyloseq::sample_data(physeq_subpop_glom)$Population)
      a <- adiv %>%
        gather(key = metric, value = value, c("Observed")) %>%
        mutate(metric = factor(metric, levels = c("Observed"))) %>%
        ggplot(aes(x = Population, y = value)) +
        geom_boxplot(outlier.color = NA) +
        geom_jitter(aes(color = Population), height = 0, width = .2) +
        labs(x = "Population", y = "Observed Richness", title = paste("Population: ", p, ", Taxa Level: ", level,sep = "")) +
        facet_wrap(~ metric, scales = "free") +
        theme(legend.position="none")
      a <- a + stat_compare_means()
      print(a)
      }

    adiv <- data.frame(
      "Observed Richness" = phyloseq::estimate_richness(physeq_subpop_glom, measures = "Observed"),
      "Quantile" = phyloseq::sample_data(physeq_subpop_glom)$Quantile)
    a <- adiv %>%
      gather(key = metric, value = value, c("Observed")) %>%
      mutate(metric = factor(metric, levels = c("Observed"))) %>%
      ggplot(aes(x = Quantile, y = value)) +
      geom_boxplot(outlier.color = NA) +
      geom_jitter(aes(color = Quantile), height = 0, width = .2) +
      labs(x = "Average Biomass Quantiles", y = "Observed Richness", title = paste("Population: ", p, ", Taxa Level: ", level,sep = "")) +
      facet_wrap(~ metric, scales = "free") +
      theme(legend.position="none")
    a <- a + stat_compare_means()
    print(a)
    
    ###################
    ###
    ### PCA
    ###
    ###################
    
    ## CLR Transform and Ordinate
    temp <- microbiome::transform(physeq_subpop_glom, 'clr')
    ord_clr <- phyloseq::ordinate(temp, "RDA", distance = "euclidean")

    # Generate distance matrix
    temp_clr_dist_matrix <- phyloseq::distance(temp, method = "euclidean")

    if (p == pops[1]){
      ## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
      #ADONIS test
      print(paste("perMANOVA by Population: ", p, ", Taxa Level: ", level,sep = ""))
      ado <- vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Population)
      print(ado)

      # Check if ADONIS confounded by differences in dispersion
      dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Population)
      print(dispr)
      print(permutest(dispr))

      #Scale axes and plot ordination
      clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
      clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
      temp2 <- phyloseq::plot_ordination(physeq_subpop_glom, ord_clr, type="samples", color="Population", title = paste("Population: ", p, ", Taxa Level: ", level, " Pr(>F): ", ado$aov.tab$`Pr(>F)`[1], sep = "")) +
        geom_point(size = .5) +
        coord_fixed(clr2 / clr1) +
        stat_ellipse(aes(group = Population), linetype = 2) +
        scale_color_manual(values=c("#78D2FF", "#18D474", "#007527", "#4222E0"))
      #dev.new()
      print(temp2)
      }

    ## Check if samples cluster beyond that expected by sampling variability w/ PERMANOVA (via vegan adonis)
    #ADONIS test
    print(paste("perMANOVA by Quantile for Population: ", p, ", Taxa Level: ", level,sep = ""))
    ado <- vegan::adonis(temp_clr_dist_matrix ~ phyloseq::sample_data(temp)$Quantile)
    print(ado)

    # Check if ADONIS confounded by differences in dispersion
    dispr <- vegan::betadisper(temp_clr_dist_matrix, phyloseq::sample_data(temp)$Quantile)
    print(dispr)
    print(permutest(dispr))

    #Scale axes and plot ordination
    clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
    clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
    temp2 <- phyloseq::plot_ordination(physeq_subpop_glom, ord_clr, type="samples", color="Quantile", title = paste("Population: ", p, ", Taxa Level: ", level, " Pr(>F): ", ado$aov.tab$`Pr(>F)`[1], sep = "")) +
      geom_point(size = 2) +
      coord_fixed(clr2 / clr1) +
      stat_ellipse(aes(group = Quantile), linetype = 2) +
      scale_color_manual(values=c("#78D2FF", "#18D474", "#007527", "#4222E0"))
    print(temp2)
    
    ###################
    ###
    ### ZETA DIVERSITY 
    ###
    ###################
    
    for (quantile in quantiles){
      temp = get(quantile)
      z <- as.data.frame(temp@otu_table)
      z[z >= 1] <- 1 #convert to presence/absence
      zeta_decline <- Zeta.decline.ex(data.spec = z, orders = 1:zeta_order, plot = FALSE)

      if(zeta_decline$aic[2,2] < zeta_decline$aic[1,2]){
        temp = "Power Law"
      } else {
        temp = "Exponential"
      }
      if(abs(zeta_decline$aic[2,2] - zeta_decline$aic[1,2]) < 2){
        temp = "Not Significant (AIC < 2)"
      } else {
        temp = temp
      }
      temp2 <- c(p, level, match(quantile,quantiles), temp, zeta_order)
      z_decline <- rbind(temp2, z_decline)
      
      for (i in 2:5){
        dec <- Zeta.ddecay(xy = x.dist.mat, data.spec = z, order = i, distance.type = "custom", dist.custom = x.dist.mat, normalize = "Jaccard", plot = FALSE)
        slope <- coef(dec$reg)[2]
        dec <- dec$reg$model
        dec_name <- paste("dec",i,sep = "")
        assign(dec_name,dec)
        temp <- c(p, level, match(quantile, quantiles), i, slope)
        z_decay <- rbind(z_decay, temp)
      }

      zp <- ggplot() +
        stat_smooth(data = dec2, method = "glm", aes(distance.reg, zeta.val.reg, color = "2")) +
        stat_smooth(data = dec3, method = "glm", aes(distance.reg, zeta.val.reg, color = "3")) +
        stat_smooth(data = dec4, method = "glm", aes(distance.reg, zeta.val.reg, color = "4")) +
        stat_smooth(data = dec5, method = "glm", aes(distance.reg, zeta.val.reg, color = "5")) +
        labs(y = "Zeta diversity", x = "Genetic Distance (Euclidean)", title = paste("Zeta Decay. Population: ",p," Quantile: ", match(quantile,quantiles), " Taxa Level: ", level, sep = "")) +
        scale_colour_manual(name = "Zeta Order", values = c("2" = "#2EFBFB", "3"="#25D8D8", "4"="#008686","5"="#005F5F"))
      print(zp)

     }
    colnames(z_decline) <- c("Population","Taxa Level","Quantile (Biomass)","Regression Model", "Zeta Order")
    colnames(z_decay) <- c("Population","Taxa Level", "Quantile (Biomass)", "Zeta Order", "Decay Slope")

    ###################
    ###
    ### SpiecEasi
    ###
    ###################
    
    if(level != "Class"){
      for(quantile in quantiles){
        temp = get(quantile)
        se.mb <- spiec.easi(temp, method='mb')#, lambda.min.ratio=1e-2, nlambda =10)#, nlambda = 100)#, nlambda = 20)#, nlambda = 10)#,lambda.min.ratio=1e-2, nlambda = 10)
        print(paste(p, level, quantile, sep = " "))
        print(getStability(se.mb))
        optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
        edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]

        optbeta <- symBeta(getOptBeta(se.mb))
        edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
        edge_weights <- edge_weights/max(edge_weights)
        ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
                             vertex.attr=list(name=taxa_names(temp)),
                             edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))

        TotalNodes = gorder(ig2.mb)
        TotalEdges = gsize(ig2.mb)
        PositiveEdges = sum(optbeta@x>0)
        NegativeEdges = sum(optbeta@x<0)
        PropPosNeg = PositiveEdges/NegativeEdges
        PosTotal = PositiveEdges/TotalEdges
        NegTotal = NegativeEdges/TotalEdges
        AvgPathLength = average.path.length(ig2.mb) #Shortest paths between vertices
        Cluster = cluster_louvain(ig2.mb) #Louvain method is an algorithm to detect communities in large networks. It maximizes a modularity score for each community, where the modularity quantifies the quality of an assignment of nodes to communities.This means evaluating how much more densely connected the nodes within a community are, compared to how connected they would be in a random network.
        Modularity = modularity(Cluster) #Strength of division of network into modules
        DegreeDist = degree(ig2.mb)
        AvgDegree = mean(DegreeDist)
        ClusterCoeff = transitivity(ig2.mb) #Transitivity ("clustering coefficient") measures the probability that the adjacent vertices of a vertex are connected.

        #degree heterogeneity
        temp2 = 0
        for(i in 1:length(table(DegreeDist))){
          temp = sum((1-(table(DegreeDist)[i]/sum(table(DegreeDist))))^2)
          temp2 = temp2+temp
        }
        HetIndex = sqrt((1/TotalNodes)*temp2)
        MaxHetInt = sqrt(1 - (3/TotalNodes) + ((TotalNodes+2)/(TotalNodes^3)))
        HetMeasure = HetIndex/MaxHetInt

        temp <- c(p, level, quantile, TotalNodes, TotalEdges, PositiveEdges, NegativeEdges, PropPosNeg, PosTotal, NegTotal, AvgPathLength, Modularity, AvgDegree, HetMeasure, ClusterCoeff)
        pos_and_neg_edges <- rbind(pos_and_neg_edges, temp)
        colnames(pos_and_neg_edges) <- c("Population", "Taxa Level", "Biomass Quantile", "Total Nodes","Total Edges", "Positive Edges", "Negative Edges", "Pos/Neg","Pos/Total","Neg/Total", "Avg Path Length", "Modularity", "Avg Degree", "Heterogeneity", "Clustering Coefficient")

        ###HUB SCORES
        if(level == "Genus"){
          temp <- as.data.frame(head(sort(hub_score(ig2.mb)$vector, decreasing = TRUE), 10))
          colnames(temp)[1] <-  "Hub Score"
          temp <- tibble::rownames_to_column(temp, "OTU")
          temp$Population = p
          temp$Quantile = quantile
          temp$Level = level

          hubscores <- rbind(hubscores, temp)

        }
        ### NETWORK GRAPH

        df <- igraph::as_data_frame(ig2.mb, 'both')
        check <- as.data.frame(taxmat[,"Phylum"])
        colnames(check) <- "Phylum"
        check$id <- row.names(check)

        df$vertices <- df$vertices %>%
          left_join(check, c('name'='id'))

        updated_g <- graph_from_data_frame(df$edges,
                                           directed = F,
                                           vertices = df$vertices)
        updated_g$layout <- layout_with_dh

        coul <- scales::viridis_pal(option = "turbo")(length(unique(tax_table(get(quantile))[,2])))

        edge_function <- paste(edge_cols,edge_weights,')',sep = "")
        edge_col_val <- c()
        for(i in 1:length(edge_function)){
          edge_col_val[i] <- eval(parse(text=edge_function[i]))
        }

        plot(updated_g, edge.color = edge_col_val, vertex.size = ((hub_score(ig2.mb)$vector)*10)+1, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Phylum")))],vertex.frame.color = "black", main = paste("Quantile:",quantile,", Pop:", p, ", Level:", level, sep = " "))

        legend("left", title = "Phylum", legend=levels(as.factor(V(updated_g)$Phylum))  , col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 1 , horiz = FALSE, inset = c(-0.15, 0))

        sizeCut<- c(0.2,0.4,0.6,0.8,1.0)
        sizeCutScale <- sizeCut*10+1
        a <- legend('right',title = "Hub Score", legend=unique(sizeCut), pt.cex= sizeCutScale, inset = c(-0, 0), bty = "n", y.intersp=1.1)
        x <- (a$text$x + a$rect$left) / 2
        y <- a$text$y
        symbols(x,y,circles=sizeCutScale/200,inches=FALSE,add=TRUE,bg='black')


      }
    }
    
  }
}

