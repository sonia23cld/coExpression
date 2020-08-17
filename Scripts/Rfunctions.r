# filter genes with too many missing data and/or zero variance
geneFilter <- function(exprData)
{
  gsg <- goodSamplesGenes(exprData, verbose = 3)
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    {printFlush(paste("Removing genes:", paste(names(expr6C)[!gsg$goodGenes], collapse = ", ")))}
    if (sum(!gsg$goodSamples)>0)
    {printFlush(paste("Removing samples:", paste(rownames(expr6C)[!gsg$goodSamples], collapse = ", ")))}
    # Remove the offending genes and samples from the data:
    return(exprData[gsg$goodSamples, gsg$goodGenes])
  }
  else
  {return(exprData)}
}


sampleTree <- function(exprData, label = sample)
{
  rownames(exprData) <- samples[match(rownames(exprData), samples$sample), label]
  sampleTree = hclust(dist(exprData), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  # clusters according accession
}

plotSoftThresholdChoices <- function(exprData)
{
  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft <- pickSoftThreshold(exprData, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

# Module Connectivity calculation. Every value is averaged based on the dimension of the module 
summarizeConnectivity <- function(expr, modules)
{
  moduleConnectivity <- data.frame('Module name'= character(), 'KWithin'= numeric(), 'KOut'= numeric(), 'Quality module'= numeric())
  for (module in modules)
  {
    geneIdx <- colnames(expr$data)[expr$mergedLabels == module]
    Kwithin <- mean(expr$geneConnectivity$kWithin[rownames(expr$geneConnectivity) %in% geneIdx])
    KOut <- mean(expr$geneConnectivity$kOut[rownames(expr$geneConnectivity) %in% geneIdx])
    Module.Quality <- Kwithin/KOut
    lineiwant <- data.frame('Module name'= module, 'KWithin'= Kwithin, 'KOut'= KOut, 'Quality module'= Module.Quality)
    moduleConnectivity <- rbind(moduleConnectivity, lineiwant)
  }
  return(moduleConnectivity)
}

#Relative module connectivity calculation. 
RelativeConnectivity <- function(expr, modules)
{
  module.relativeConnectivity <- data.frame('Module name'= character(), 'KWithin'= numeric(), 'KOut'= numeric(), 'KDiff'=numeric(), 'Quality module'= numeric())
  for (module in modules)
  {
    geneIdx <- colnames(expr$data)[expr$mergedLabels == module]
    Kwithin<- c()
    KOut<- c()
    for (gene in geneIdx) {
      Kwithin <- c(Kwithin, expr$geneConnectivity$kWithin[rownames(expr$geneConnectivity) == gene])
      KOut <- c(KOut, expr$geneConnectivity$kOut[rownames(expr$geneConnectivity) == gene])
    }
    relative.Kwithin <- Kwithin/(length(geneIdx)-1)
    relative.KOut <- KOut/(length(colnames(expr$data))-(length(geneIdx)))
    KWithin.module<- mean(relative.Kwithin)
    KOut.module <- mean(relative.KOut)
    KDiff.module<- KWithin.module - KOut.module
    Module.Quality <- KWithin.module/KOut.module
    lineiwant <- data.frame('Module name'= module, 'KWithin'= KWithin.module, 'KOut'= KOut.module, 'KDiff'= KDiff.module, 'Quality module'= Module.Quality)
    module.relativeConnectivity <- rbind(module.relativeConnectivity, lineiwant)
  }
  return(module.relativeConnectivity)
}

#Centrality. (The result is the same you would have with the Relative KOut)
Centrality<- function(all.expr, modules) 
{
  Centrality<- c()
  for (module in all.expr$moduleConnectivity$Module.name)
  {
    mean_KOut<- c()
    geneIdx <- colnames(all.expr$data)[all.expr$mergedLabels == module]
    n_genes_out<- as.numeric(ncol(all.expr$data)-length(geneIdx))
    for (gene in geneIdx) {
      KOut <- all.expr$geneConnectivity$kOut[rownames(all.expr$geneConnectivity) %in% gene]
      mean_KOut<- c(mean_KOut, KOut/n_genes_out)
    }
    Centrality<- c(Centrality, mean(mean_KOut))
  }
  return(Centrality)
}

#create color-coded table of the intersection counts
# Truncate p values smaller than 10^{-50} to 10^{-50}
Heatmap_overlap<- function(pTable, CountTbl, expr) {
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
  pTable[pTable>50 ] = 50 
  # Marginal counts (really module sizes)
  ModTotal<- list(apply(CountTbl, 1, sum), apply(CountTbl, 2, sum))
  # Actual plotting
  sizeGrWindow(10,7)
  #pdf(file = "/Volumes/nordborg/pub/forPieter/WGCNA/Results/Overlap modules 6 vs16.pdf", wi = 10, he = 7);
  par(mfrow=c(1,1));
  par(cex = 1.0);
  par(mar=c(8, 10.4, 2.7, 1)+0.3)
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  labeledHeatmap(Matrix = pTable,
                 xLabels = paste(" ", expr[[2]]$Modules),
                 yLabels = paste(" ", expr[[1]]$Modules),
                 colorLabels = TRUE,
                 xSymbols = paste("16 ", expr[[2]]$Modules, ": ", ModTotal[[2]], sep=""),
                 ySymbols = paste("6 ", expr[[1]]$Modules, ": ", ModTotal[[1]], sep=""),
                 textMatrix = CountTbl,
                 colors = greenWhiteRed(100)[50:100],
                 main = "Correspondence of 6 set-specific and 16 set-specific modules",
                 cex.text = 0.3, cex.lab = 0.3, setStdMargins = FALSE);
}

Table.overlap<- function(expr, pTable, Threshold, CountTbl) {
  for (mod6 in 1:length(expr[[1]]$Modules)) {
    for (mod16 in 1:length(expr[[2]]$Modules)) {
      if (pTable[mod6, mod16] > Threshold) {
        lineiwant<- data.frame('Name module 16C'= expr[[2]]$Modules[mod16], 
                               'Name module 6C'= expr[[1]]$Modules[mod6], 
                               'N_totgenes16'= count(expr[[2]]$mergedColors == expr[[2]]$Modules[mod16]), 
                               'N_totgenes6'= count(expr[[1]]$mergedColors == expr[[1]]$Modules[mod6]), 
                               'N_overlapped_genes'= CountTbl[mod6, mod16],
                               'p-value'=pTable[mod6,mod16])
        ClearTable<- rbind(ClearTable, lineiwant) 
        next }
    }
  }
  return(ClearTable)
}

#Function Go enrichment
Top.Go<- function (genesinmodule, go_ids, genesuniverse, gene_2_GO) {
  # remove any candidate genes without GO annotation
  keep = genesinmodule %in% go_ids[,2]
  keep =which(keep==TRUE)
  mygenes=genesinmodule[keep]
  
  # make named factor showing which genes are of interest
  geneList=factor(as.integer(genesuniverse %in% mygenes))
  names(geneList)= genesuniverse
  #make topGO data object
  GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
  
  # define test using the classic algorithm with fisher 
  classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
  # define test using the weight01 algorithm (default) with fisher
  weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
  
  # generate a table of results: GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(GOdata)
  all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))
  #performing BH correction on our p values
  p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
  
  # create the file with all the statistics from GO analysis
  all_res_final=cbind(all_res,p.adj)
  all_res_final=all_res_final[order(all_res_final$p.adj),]
  
  Results<- as.data.frame(all_res_final[1:50,])
  return(Results)
}

Genes.behaviour<- function(xrange, yrange, Tableplot, expr) {
  #set up the plot
  plot(NA, NA, xlim=xrange, ylim = yrange, xaxt = 'n', type='n', xlab='Name sample', ylab='Expression level', main=paste('Expression level genes in module', module,  names(expr), sep = ' '))
  axis(1, at=1:12, labels = rownames(Tableplot), cex.axis=0.5)
  colors <- rainbow(ncol(Tableplot))
  #legend('topright', inset = c(-0.2,-0.1), colnames(Tableplot), fill =colors, cex=0.5)
  #add lines
  for (i in 1:ncol(Tableplot)) {
    lines(Tableplot[i], type='l', lwd=1.5, col=colors[[i]])
  }
}

#find snps in promoter region
gene_SNPmatrix_promoter<- function(chr_regions, Chr, coordinates_gene){
  indexes<- chr_regions[, Chr]
  positions<- h5read('/Volumes/nordborg/pub/forPieter/WGCNA/all_chromosomes_binary_Rcompatible.hdf5', '/positions', index = list((indexes[1]+1):indexes[2]))
  gene_indexes_before<- which(positions %in% c((coordinates_gene[[1]]-10000):(coordinates_gene[[1]]-1)))
  gene_indexes_before<- gene_indexes_before + indexes[1]
  gene_indexes_after<- which(positions %in% c((coordinates_gene[[2]]+1):(coordinates_gene[[2]]+10000)))
  gene_indexes_after<- gene_indexes_after + indexes[1]
  SNPsmatrix_before<- h5read('/Volumes/nordborg/pub/forPieter/WGCNA/all_chromosomes_binary_Rcompatible.hdf5', '/snps', index = list(1:1135, min(gene_indexes_before):max(gene_indexes_before)))
  SNPsmatrix_after<- h5read('/Volumes/nordborg/pub/forPieter/WGCNA/all_chromosomes_binary_Rcompatible.hdf5', '/snps', index = list(1:1135, min(gene_indexes_after):max(gene_indexes_after)))
  SNPsmatrix<- cbind(SNPsmatrix_before, SNPsmatrix_after)
  rownames(SNPsmatrix)<- h5read('/Volumes/nordborg/pub/forPieter/WGCNA/all_chromosomes_binary_Rcompatible.hdf5', '/accessions')
  colnames(SNPsmatrix)<- c(positions[positions %in% c((coordinates_gene[[1]]-10000):(coordinates_gene[[1]]-1))], positions[positions %in% c((coordinates_gene[[2]]+1):(coordinates_gene[[2]]+10000))])
  return(SNPsmatrix)}                                                                                                           

#minor allele frequency cut-off
MAF<- function(SNPsmatrix) {
  n0<- apply(SNPsmatrix, 2, function(x) sum(x == 0))
  allele_freq<- unlist(lapply(n0, function(x) {(x * 100)/nrow(SNPsmatrix)}))
  cut_off<- c(which(allele_freq <10), which(allele_freq >90))
  SNPsmatrix<- as.data.table(SNPsmatrix)[, (cut_off) := NULL]}


