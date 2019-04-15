# WGCNA


options(stringsAsFactors = FALSE)


library(tximport)
library(DESeq2)
source('Scripts/Rfunctions.r')

setwd('/Volumes/nordborg/user/pieter.clauw/Documents/Experiments/UltimateQandD/')
# DATA 
samples<- read.table('/Volumes/nordborg/pub/forPieter/WGCNA/samples.txt', header = T, comment.char = '', sep = '\t')
Araport11 <- read.table('/Volumes/nordborg/pub/forPieter/WGCNA/Araport11_GFF3_genes_transposons.201606.ChrM_ChrC_FullName.gtf')
# FUNCTION



# PREPARATION
samples <- samples[!samples$temperature %in% c('6C?', '16C?'), ]
samples <- samples[samples$genoCheck == 'ok', ]
samples$accession <- as.factor(samples$accession)
samples$temperature <- as.factor(samples$temperature)
samples$replicate <- as.factor(samples$replicate)
files <- file.path('/Volumes/nordborg/pub/forPieter/WGCNA/SalmonQuantification/', samples$basename, '' ,'_quasiMap_pseudoG/quant.sf', fsep = '')
names(files) <- samples$sample

colnames(Araport11) <- c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
Araport11$attributes <- as.character(Araport11$attributes)

tx2gene <- data.frame('transcriptID' = Araport11[,10], 'geneID' = Araport11[,13])

txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, dropInfReps = T)
dds.full <- DESeqDataSetFromTximport(txi, colData = samples, design = ~  replicate + accession + temperature + accession:temperature)
dds <- estimateSizeFactors(dds.full)
idx <- rowSums(counts(dds) >= 10 ) >= 19.80
dds<- dds[idx,]

# variance stabilisation transform of count data
dds.varStab <- varianceStabilizingTransformation(dds, blind = F)

#remove batch effect
library(limma)
assay(dds.varStab) <- removeBatchEffect(assay(dds.varStab), dds.varStab$replicate)

counts.varStab <- t(assay(dds.varStab))

detach("package:DESeq2", unload = T)
detach("package:tximport", unload = T)
library(WGCNA)
# WGCNA
expr6C <- counts.varStab[as.character(samples$sample[samples$temperature == '6C']), ]
expr16C <- counts.varStab[as.character(samples$sample[samples$temperature == '16C']), ]

#check quantile scatterplots to make sure there are no systematic shifts between samples
multiExpr<-list('expr16C'=expr16C, 'expr6C'=expr6C)
for (i in c(1, 2)) {
  readmatrix<- multiExpr[[i]]
  for (j in 1:nrow(readmatrix)) {
    sampleiwant<- t(readmatrix[j,])
    pdf(file = paste('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Check normalization/', row.names(readmatrix)[j], '.pdf', sep = ''), wi = 14, he = 8)
    qqnorm(sampleiwant, pch = 1, frame = FALSE, main = paste(names(multiExpr[i]), 'sample', row.names(readmatrix)[j], sep = ' '))
    qqline(sampleiwant, col = "red", lwd = 2)
    dev.off()
  }
}


# filter genes with too many missing values and zero variance
expr6C.filt <- geneFilter(expr6C)
expr16C.filt <- geneFilter(expr16C)

# check clustering according accessions
sampleTree(expr6C.filt, label = 'accession')
sampleTree(expr16C.filt, label = 'accession')
sampleTree(expr16C.filt, label = 'replicate')
# acessions much more clearly clustered at 6C compared to 16C


# decide on softThreshold
plotSoftThresholdChoices(expr6C.filt)
# pick softhtThreshold 12 for 6C data
plotSoftThresholdChoices(expr16C.filt)
# pick softThreshold 10

# clustering
expr6C.net <- blockwiseModules(expr6C.filt, maxBlockSize = 20000,
                         power = 12, TOMType = "unsigned", minModuleSize = 10,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = 'expr6c.TOM.blockwise',
                         verbose = 3)

expr16C.net <- blockwiseModules(expr6C.filt, maxBlockSize = 20000,
                                power = 10, TOMType = "unsigned", minModuleSize = 10,
                                reassignThreshold = 0, mergeCutHeight = 0.25,
                                numericLabels = TRUE,
                                saveTOMs = TRUE,
                                saveTOMFileBase = 'expr6c.TOM.blockwise',
                                verbose = 3)





# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(expr6C.net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(expr6C.net$dendrograms[[1]], mergedColors[expr6C.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


GOenr = GOenrichmentAnalysis(mergedColors, allLLIDs, organism = "mouse", nBestP = 10)


# export results
for (module in unique(expr16C.net$colors))
{
  # Select module probes
  modGenePresence <- (expr16C.net$colors==module)
  # Get their entrez ID codes
  modGenes = colnames(expr16C.filt)[modGenePresence]
  # Write them into a file
  fileName = paste("Results/Transcriptome/WGCNA/16C/wgcna_16C_GenesinModule_", module, ".txt", sep="");
  write.table(as.data.frame(modGenes), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

