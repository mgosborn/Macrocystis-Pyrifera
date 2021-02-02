library(vegan)
library(plotly)

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

veg_sampledata <- pssd2veg(physeq_merged_transformed)
veg_otu <- psotu2veg(physeq_merged_transformed)

# Calculating relative abundance and creating new dataframe with relative abundance data
veg_otu.rel <- decostand(veg_otu, method = "total")

# Calculate distance matrix
veg_otu_distmat <- vegdist(veg_otu.rel, method = "bray")

# Creating easy to view matrix and writing .csv
veg_otu_distmat <- as.matrix(veg_otu_distmat, labels = T)
write.csv(veg_otu_distmat, "veg_otu_distmat.csv")

# Running NMDS in vegan (metaMDS)
veg_otu_NMS <-
  metaMDS(veg_otu_distmat,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

OTUxyz = scores(veg_otu_NMS, display="sites")
plot_ly(x=OTUxyz[,1], y=OTUxyz[,2], z=OTUxyz[,3], text= "Bray Curtis NMDS", type="scatter3d", mode="markers", color=veg_sampledata$Population)
plot_ly(x=OTUxyz[,1], y=OTUxyz[,2], z=OTUxyz[,3], text= "Bray Curtis NMDS", type="scatter3d", mode="markers", color = veg_sampledata$Stipe_Rank)

dev.new()
par(mfrow=c(1,2))
#Axis 1 and 2 (x and y)
plot(OTUxyz[,1], OTUxyz[,2], main="Bray-Curtis 1:2", pch=20, col=c("darkslategray3","darksalmon","lightsteelblue3","mediumorchid1")[veg_sampledata$Population])
legend(-.37, -.36, legend=c("AQ","CB","CI","LC"), col=c("darkslategray3","darksalmon","lightsteelblue3","mediumorchid1"), pch=20)
#Axis 1 and 3 (x and z)
plot(OTUxyz[,1], OTUxyz[,3], main="Bray-Curtis 1:3", pch=20, col=c("darkslategray3","darksalmon","lightsteelblue3","mediumorchid1")[veg_sampledata$Population])


# Shepards test/goodness of fit
goodness(veg_otu_NMS) # Produces a results of test statistics for goodness of fit for each point
stressplot(veg_otu_NMS) # Produces a Shepards diagram

# Plotting points in ordination space
plot(veg_otu_NMS, "sites")   # Produces distance 
orditorp(veg_otu_NMS, "sites")   # Gives points labels
