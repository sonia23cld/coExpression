setwd('~/Network/')
options(stringsAsFactors = FALSE)
library(WGCNA)
load("data/multiExpr_newparameters.RData")
exprSize = checkSets(multiExpr)
nGenes = exprSize$nGenes

softPower = 7

#for 16C
adjacency16C = abs(bicor(multiExpr[[1]]$data, use = "p", maxPOutliers = 0.10))^softPower
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency16C)
dissTOM = 1-TOM 
# Call the hierarchical clustering function
geneTree16C = hclust(as.dist(dissTOM), method = "average")

# set the minimum module size 
minModuleSize = 20
# Module identification using dynamic tree cut:
dynamicMods16C = cutreeDynamic(dendro = geneTree16C, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)

# Convert numeric lables into colors
dynamicColors16C = labels2colors(dynamicMods16C)

softPower=7
#for 6C
adjacency6C = abs(bicor(multiExpr[[2]]$data, use = "p", maxPOutliers = 0.10))^softPower
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency6C)
dissTOM = 1-TOM 
# Call the hierarchical clustering function
geneTree6C = hclust(as.dist(dissTOM), method = "average")

# set the minimum module size 
minModuleSize = 20
# Module identification using dynamic tree cut:
dynamicMods6C = cutreeDynamic(dendro = geneTree6C, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)

# Convert numeric lables into colors
dynamicColors6C = labels2colors(dynamicMods6C)

save(dynamicMods16C, dynamicColors16C, geneTree16C, dynamicMods6C, dynamicColors6C, geneTree6C, file = '/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft40/Set-specific network like consensus40.RData')
