#for consensus
sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#merging modules
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "", cex=0.5)
abline(h=0.15, col = "red")

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.15, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

save(moduleColors, moduleLabels, consMEs, file = '/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft10/Consensus_aftermerging10_015.RData')


#for specific
geneTree<-geneTree6C
dynamicColors<- dynamicColors6C
dynamicMods<- dynamicMods6C
multiExpression<- multiExpr[[2]]$data
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(multiExpression, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "", cex=0.5)

MEDissThres = 0.15
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(multiExpression, dynamicMods, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
moduleLabels = merge$colors;
mergedColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# Rename to moduleColors
moduleColors = mergedColors

MEs = mergedMEs

moduleColors6C<- moduleColors
moduleLabels6C<- moduleLabels
MEs6C<- MEs

save(moduleColors16C, moduleLabels16C, MEs16C, moduleColors6C, moduleLabels6C, MEs6C, file = '/Volumes/nordborg/pub/forPieter/WGCNA/Results/New data/Setspecific_aftermerging07_015.RData')


#check connectivty

#method from expression data
connectivity<-intramodularConnectivity.fromExpr(multiExpr[[1]]$data, moduleLabels16C, 
                                  corFnc = "bicor", corOptions = ("use = 'p', maxPOutliers = 0.10"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'average'",
                                  networkType = "unsigned", power = 7,
                                  scaleByMax = FALSE,
                                  ignoreColors = if (is.numeric(moduleLabels16C)) 0 else "grey",
                                  getWholeNetworkConnectivity = TRUE)


#method from adjacency matrix
connectivity<-intramodularConnectivity(adjacency16C, moduleLabels16C, scaleByMax = FALSE)



#save(connectivity16C_07, connectivity6C_07, file = '/Volumes/nordborg/pub/forPieter/WGCNA/Results/New data/Setspecific_connectivity07.RData')


#connectivity for modules
connectivity<- connectivity6C_07
modulecolor<-moduleColors6C
multiExpression<- multiExpr[[2]]$data
connectivity[is.na(connectivity)] <- 0 
connectivity$gene<- (colnames(multiExpression))
connectivitymodules<- data.frame(matrix(0, ncol = 5, nrow = length(unique(modulecolor))))
colnames(connectivitymodules)<- colnames(connectivity)
colnames(connectivitymodules)[which(names(connectivity) == 'gene')]<- c('Name_module')
connectivitymodules$Quality_check<- NA
connectivitymodules<- connectivitymodules[, -c(1,4)]
connectivitymodules$Name_module<- unique(modulecolor)
for (i in 1:nrow(connectivitymodules)) {
   namemodule<- connectivitymodules$Name_module[[i]]
   genesinmodule<- colnames(multiExpression)[modulecolor == namemodule]
   connectivitymodules[connectivitymodules$Name_module %in% namemodule, 'kOut'] <- (sum(connectivity[connectivity$gene %in% genesinmodule,]$kOut))/length(genesinmodule)
   connectivitymodules[connectivitymodules$Name_module %in% namemodule, 'kWithin'] <- (sum(connectivity[connectivity$gene %in% genesinmodule,]$kWithin))/length(genesinmodule)
   connectivitymodules[connectivitymodules$Name_module %in% namemodule, 'Quality_check']<- (connectivitymodules[connectivitymodules$Name_module %in% namemodule,]$kWithin)/(connectivitymodules[connectivitymodules$Name_module %in% namemodule,]$kOut)
  }

connectivitymodules[is.na(connectivitymodules)] <- 0 
module_qualitycheck<- sum(connectivitymodules$Quality_check)/nrow(connectivitymodules)

