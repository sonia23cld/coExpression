#Network construction and consensus module detection
##Use data prepared in WGCNA script
library(WGCNA)
setwd('/Volumes/nordborg/pub/forPieter/WGCNA/')
nSets = 2
setLabels = c("16C", "6C")

#Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(expr16C))
multiExpr[[2]] = list(data = as.data.frame(expr6C))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK #it's false so we have to correct it



if (!gsg$allOK) {
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
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



#cluster on distance every sample
sampleTrees = list()
for (set in 1:nSets) {
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
#plot dendogram
#pdf(file = "SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets) {
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)}
dev.off()
collectGarbage()


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets) {
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2)[[2]])
}
collectGarbage()
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets) {
  for (col in 1:length(plotCols)) {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets) {
  if (set==1) {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)", ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1) {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers, cex=cex1, col=colors[set])
  } else {
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set])
  }
  if (col==1) {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) 
  } else {
    legend("topright", legend = setLabels, col = colors, pch = 20) 
  }
}
dev.off()

#network construction
bnet = blockwiseConsensusModules(multiExpr, maxBlockSize = 20000,
                       power = 14, TOMType = "unsigned", minModuleSize = 30, deepSplit = 2, pamRespectsDendro = FALSE,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       verbose = 5)


#see the result
consMEs = bnet$multiMEs;
moduleLabels = bnet$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = bnet$dendrograms[[1]];
bwLabels = matchLabels(bnet$colors, moduleLabels, pThreshold = 1e-7);
bwColors = labels2colors(bwLabels)

table(bwLabels)
bwLabels

#plot the dendrogram and module color
# Here we show a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
#pdf(file = "Plots/BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(bnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks) {
  plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
}

dev.off()


