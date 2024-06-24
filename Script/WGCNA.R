if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("preprocessCore", "impute" ))

library(dplyr)
library(GEOquery)
library(tidyverse) 
library(DESeq2)
library(data.table)
library(WGCNA)
library(CorLevelPlot)
library(ggplot2)
library(gridExtra) 
library(dynamicTreeCut)

options(stringsAsFactors = FALSE)

## tranpose data - swap genes to columns and samples to rows
trans_111006 <- t(vsd_111006_df)
trans_111010 <- t(vsd_111010_df)
trans_111016 <- t(vsd_111016_df)
trans_113165 <- t(vsd_113165_df)
trans_120642 <- t(vsd_120642_df)
trans_126865 <- t(vsd_126865_df)
trans_144304 <- t(vsd_144304_df)
trans_151066 <- t(vsd_151066_df)
trans_159217 <- t(vsd_159217_df)
trans_163434 <- t(vsd_163434_df)
trans_164471 <- t(vsd_164471_df)
trans_165630 <- t(vsd_165630_df)
trans_174106 <- t(vsd_174106_df)
trans_175495 <- t(vsd_175495_df)
trans_226151 <- t(vsd_226151_df)
trans_235781 <- t(vsd_235781_df)
trans_242202 <- t(vsd_242202_df)

## choose a set of threshold powers 
power_v <- c(seq(1,10, by=1), seq(12,30, by=2)) 

## soft threshold 
##----call the network topology analysis function-111006 determine soft threshold-----
soft_111006 <- pickSoftThreshold(trans_111006, powerVector = power_v,
                                networkType = "signed", verbose = 2)
soft_111006_data <- soft_111006$fitIndices

## plot the results to visualize the data
a1<- ggplot(soft_111006_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a2<- ggplot(soft_111006_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a1, a2, nrow = 2)
## noting my soft threshold 
softthresh_111006 = 8



##-----soft threshold 111010----
soft_111010 <- pickSoftThreshold(trans_111010, powerVector = power_v,
                                 networkType = "signed", verbose = 2)
soft_111010_data <- soft_111010$fitIndices

## plot the results to visualize the data
a3<- ggplot(soft_111010_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a4<- ggplot(soft_111010_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a3, a4, nrow = 2)
softthresh_111010 = 9



##--------soft threshold 111016----

soft_111016 <- pickSoftThreshold(trans_111016, powerVector = power_v,
                                 networkType = "signed", verbose = 2)
soft_111016_data <- soft_111016$fitIndices

## plot the results to visualize the data
a5<- ggplot(soft_111016_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a6<- ggplot(soft_111016_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity')

## combines both graphs in one box
grid.arrange(a5, a6, nrow = 2)
softthresh_111016 = 12 
 

##---- soft threshold 113165------
soft_113165 <- pickSoftThreshold(trans_113165, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_113165_data <- soft_113165$fitIndices

## plot the results to visualize the data
a7<- ggplot(soft_113165_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a8<- ggplot(soft_113165_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a7, a8, nrow = 2)
softthresh_113165 = 14

##--- soft threshold 120642------ 
soft_120642 <- pickSoftThreshold(trans_120642, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_120642_data <- soft_120642$fitIndices

## plot the results to visualize the data
a9<- ggplot(soft_120642_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a10<- ggplot(soft_120642_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a9, a10, nrow = 2)
softthresh_120642 = 10

##----soft threshold 126865------
soft_126865 <- pickSoftThreshold(trans_126865, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_126865_data <- soft_126865$fitIndices

## plot the results to visualize the data
a11<- ggplot(soft_126865_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a12<- ggplot(soft_126865_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a11, a12, nrow = 2)
softthresh_126865 = 16

##---- soft threshold 144304------
soft_144304 <- pickSoftThreshold(trans_144304, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_144304_data <- soft_144304$fitIndices

## plot the results to visualize the data
a13<- ggplot(soft_144304_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a14<- ggplot(soft_144304_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a13, a14, nrow = 2) 
softthresh_144304 = 14

##---- soft threshold 151066----
soft_151066 <- pickSoftThreshold(trans_151066, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_151066_data <- soft_151066$fitIndices

## plot the results to visualize the data
a15<- ggplot(soft_151066_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a16<- ggplot(soft_151066_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a15, a16, nrow = 2)
softthresh_151066 = 14

##--- soft threshold 159217----- 
soft_159217 <- pickSoftThreshold(trans_159217, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_159217_data <- soft_159217$fitIndices

## plot the results to visualize the data
a17<- ggplot(soft_159217_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a18<- ggplot(soft_159217_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a17, a18, nrow = 2)
softthresh_159217 = 14

##--- soft threshold 163434 ???-----
soft_163434 <- pickSoftThreshold(trans_163434, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_163434_data <- soft_163434$fitIndices

## plot the results to visualize the data
a19<- ggplot(soft_163434_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a20<- ggplot(soft_163434_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a19, a20, nrow = 2) 
softthresh_163434

##--- s0ft threshold 164471 ---- 
soft_164471 <- pickSoftThreshold(trans_164471, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_164471_data <- soft_164471$fitIndices

## plot the results to visualize the data
a21<- ggplot(soft_164471_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a22<- ggplot(soft_164471_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a21, a22, nrow = 2) 
softthresh_164471 = 14

##--- soft threshold 165630-----
soft_165630 <- pickSoftThreshold(trans_165630, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_165630_data <- soft_165630$fitIndices

## plot the results to visualize the data
a23<- ggplot(soft_165630_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a24<- ggplot(soft_165630_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a23, a24, nrow = 2) 
softthresh_165630 = 16

##--- soft threshold 174106-----
soft_174106 <- pickSoftThreshold(trans_174106, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_174106_data <- soft_174106$fitIndices

## plot the results to visualize the data
a25<- ggplot(soft_174106_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a26<- ggplot(soft_174106_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a25, a26, nrow = 2) 
softthresh_174106 = 10


##---soft threshold 175495----
soft_175495 <- pickSoftThreshold(trans_175495, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_175495_data <- soft_175495$fitIndices

## plot the results to visualize the data
a27<- ggplot(soft_175495_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a28<- ggplot(soft_175495_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a27, a28, nrow = 2) 
softthresh_175495 = 12

##--- soft threshold 226151------
soft_226151 <- pickSoftThreshold(trans_226151, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_226151_data <- soft_226151$fitIndices

## plot the results to visualize the data
a29<- ggplot(soft_226151_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a30<- ggplot(soft_226151_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a29, a30, nrow = 2) 
softthresh_226151 = 14

##---- soft threshold 235781-----
soft_235781 <- pickSoftThreshold(trans_235781, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_235781_data <- soft_235781$fitIndices

## plot the results to visualize the data
a31<- ggplot(soft_235781_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a32<- ggplot(soft_235781_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a31, a32, nrow = 2) 
softthresh_235781 = 14


##--- soft threshold 242202----
soft_242202 <- pickSoftThreshold(trans_242202, powerVector = power_v,
                                 networkType = "signed", verbose = 2) 
soft_242202_data <- soft_242202$fitIndices

## plot the results to visualize the data
a33<- ggplot(soft_242202_data, aes(Power, SFT.R.sq, label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'blue') +
  labs(x = 'Soft Threshold (Power)',
       y = 'Scale free topology model fit, signed R^2') +
  theme_grey() 

a34<- ggplot(soft_242202_data, aes(Power, mean.k., label= Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Soft Threshold (Power)', y = 'Mean Connectivity') 

## combines both graphs in one box
grid.arrange(a33, a34, nrow = 2) 
softthresh_242202 = 12




##-- co-expression similarities and adjacencies-----
adj_111006 <- adjacency(trans_111006, type = "signed", power = softthresh_111006)
adj_111010 <- adjacency(trans_111010, type = "signed", power = softthresh_111010)
adj_111016 <- adjacency(trans_111016, type = "signed", power = softthresh_111016)
adj_113165 <- adjacency(trans_113165, type = "signed", power = softthresh_113165) 
adj_120642 <- adjacency(trans_120642, type = "signed", power = softthresh_120642)
adj_126865 <- adjacency(trans_126865, type = "signed", power = softthresh_126865) 
adj_144304 <- adjacency(trans_144304, type = "signed", power = softthresh_144304)
adj_151066 <- adjacency(trans_151066, type = "signed", power = softthresh_151066)
adj_159217 <- adjacency(trans_159217, type = "signed", power = softthresh_159217) 
adj_164471 <- adjacency(trans_164471, type = "signed", power = softthresh_164471) 
adj_165630 <- adjacency(trans_165630, type = "signed", power = softthresh_165630) 
adj_174106 <- adjacency(trans_174106, type = "signed", power = softthresh_174106) 
adj_175495 <- adjacency(trans_175495, type = "signed", power = softthresh_175495) 
adj_226151 <- adjacency(trans_226151, type = "signed", power = softthresh_226151) 
adj_235781 <- adjacency(trans_235781, type = "signed", power = softthresh_235781)
adj_242202 <- adjacency(trans_242202, type = "signed", power = softthresh_242202)



##--- Topological Overlap Matrix (TOM)-----
## Turn adjacency into topological overlap 111006
TOM_111006 <- TOMsimilarity(adj_111006, TOMType = "signed") 
dissTOM_111006 <- 1 - TOM_111006
## TOM111010
TOM_111010 <- TOMsimilarity(adj_111010, TOMType = "signed")
dissTOM_111010 <- 1 - TOM_111010 
## TOM111016
TOM_111016 <- TOMsimilarity(adj_111016, TOMType = "signed")
dissTOM_111016 <- 1 - TOM_111016
## TOM113165
TOM_113165 <- TOMsimilarity(adj_113165, TOMType = "signed")
dissTOM_113165 <- 1 - TOM_113165
## TOM120642
TOM_120642 <- TOMsimilarity(adj_120642, TOMType = "signed")
dissTOM_120642 <- 1 - TOM_120642
## TOM126865
TOM_126865 <- TOMsimilarity(adj_126865, TOMType = "signed")
dissTOM_126865 <- 1 - TOM_126865
## TOM144304
TOM_144304 <- TOMsimilarity(adj_144304, TOMType = "signed")
dissTOM_144304 <- 1 - TOM_144304
## TOM151066
TOM_151066 <- TOMsimilarity(adj_151066, TOMType = "signed")
dissTOM_151066 <- 1 - TOM_151066
## TOM159217
TOM_159217 <- TOMsimilarity(adj_159217, TOMType = "signed")
dissTOM_159217 <- 1 - TOM_159217
## TOM164471
TOM_164471 <- TOMsimilarity(adj_164471, TOMType = "signed")
dissTOM_164471 <- 1 - TOM_164471
## TOM165630
TOM_165630 <- TOMsimilarity(adj_165630, TOMType = "signed")
dissTOM_165630 <- 1 - TOM_165630
## TOM174106
TOM_174106 <- TOMsimilarity(adj_174106, TOMType = "signed")
dissTOM_174106 <- 1 - TOM_174106
## TOM175495
TOM_175495 <- TOMsimilarity(adj_175495, TOMType = "signed")
dissTOM_175495 <- 1 - TOM_175495
## TOM226151
TOM_226151 <- TOMsimilarity(adj_226151, TOMType = "signed")
dissTOM_226151 <- 1 - TOM_226151
## TOM235781
TOM_235781 <- TOMsimilarity(adj_235781, TOMType = "signed")
dissTOM_235781 <- 1 - TOM_235781
## TOM242202
TOM_242202 <- TOMsimilarity(adj_242202, TOMType = "signed")
dissTOM_242202 <- 1 - TOM_242202












## clustering using TOMs ------
## call the hierarchical clustering function
geneTree_111006 <- hclust(as.dist(dissTOM_111006), method = "average")
# Plot the resulting clustering tree (dendrogram)
#create a graphic window
sizeGrWindow(12,9) 
plot(geneTree_111006, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# We like large modules, so we set the minimum module size relatively high: (use for all)
# minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods_111006 = cutreeDynamic(dendro = geneTree_111006, distM = dissTOM_111006,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)
table(dynamicMods_111006)

# Convert numeric lables into colors
dynamicColors_111006 = labels2colors(dynamicMods_111006)
table(dynamicColors_111006) 
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_111006, dynamicColors_111006, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors") 
# merge modules whose expression profiles are similar
# Calculate eigengenes
MEList_111006 = moduleEigengenes(trans_111006, colors = dynamicColors_111006)  
MEs_111006 = MEList_111006$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss_111006 = 1-cor(MEs_111006) 
# Cluster module eigengenes
METree_111006 = hclust(as.dist(MEDiss_111006), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree_111006, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres_111006 = 0.3
# Plot the cut line into the dendrogram
abline(h=MEDissThres_111006, col = "green")
# Call an automatic merging function
merge_111006 = mergeCloseModules(trans_111006, dynamicColors_111006,
                                 cutHeight = MEDissThres_111006, verbose = 3)
# The merged module colors
mergedColors_111006 = merge_111006$colors
# Eigengenes of the new merged modules:
mergedMEs_111006 = merge_111006$newMEs

## observe effect of merge on module colors - plot dendogram
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree_111006, cbind(dynamicColors_111006, mergedColors_111006),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



