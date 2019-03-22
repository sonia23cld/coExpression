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
