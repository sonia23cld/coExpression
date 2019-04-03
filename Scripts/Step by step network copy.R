#Step by step network construction


library(WGCNA)
options(stringsAsFactors = FALSE)
setLabels = c("16C", "6C")

#data input
expr6C<- read.table('/Volumes/nordborg/pub/forPieter/WGCNA/Results/data/expr6C.txt', header = TRUE, sep = '')
expr16C<- read.table('/Volumes/nordborg/pub/forPieter/WGCNA/Results/data/expr16C.txt', header = TRUE, sep = '')
nSets<- 2
multiExpr<- vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(expr16C))
multiExpr[[2]] = list(data = as.data.frame(expr6C))

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK #it's false so we have to correct it


if (!gsg$allOK) {
  # Print information about the removed genes:
  #if (sum(!gsg$goodGenes) > 0) 
    #printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
  for (set in 1:exprSize$nSets) {
    if (sum(!gsg$goodSamples[[set]])) 
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

nGenes = exprSize$nGenes
nSamples = exprSize$nSamples
nSets = checkSets(multiExpr)$nSets


#data are ready to start with analysis 
#I jump the step to find a set of soft-thresholding powers


#calculation of network adjacencies
softPower = 14
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate adjacencies in each individual data set
for (set in 1:nSets) {
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
}




# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set
for (set in 1:nSets) {
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
}


#Scaling of Topological Overlap Matrices to make them comparable across sets
# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples) 
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets) {
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1) {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
    TOM[set, ,] = TOM[set, ,]^scalePowers[set]
  }
}

#plot the 2 datasets topological overlaps before and after scaling
# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list()
for (set in 1:nSets) {
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
}
# Open a suitably sized graphics window
sizeGrWindow(6,6)
pdf(file = "~/CoExpression/TOMScaling-QQPlot.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()


#Calculation of consensus Topological Overlap
consensusTOM = pmin(TOM[1, , ], TOM[2, , ])
#Clustering and module identification
consTree = hclust(as.dist(1-consensusTOM), method = "average")
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)


#Merging of modules whose expression profiles are very similar
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs)
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average")

#merge modules under threshold
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)


#take informations you need
# Numeric module labels
moduleLabels = merge$colors
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs

#plot the gene dendrogram with both the unmerged and the merged module colors
sizeGrWindow(9,6)
pdf(file = "~/CoExpression/TOMScaling-QQPlot.pdf", wi = 6, he = 6);
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

