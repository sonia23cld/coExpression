# differential expression functions


# QC plot
QCplot <- function(QC, parameter)
{
  # colour and symbol definitions
  temperatures <- unique(QC$temperature)
  replicates <- unique(QC$replicate)
  tempCol <- data.frame('temp' = temperatures, 'col' = c('red', 'blue'))
  repSymb <- data.frame('rep' = replicates, 'symb' = c(15,16,17))
  
  QC <- QC[, c('sample', 'accession', 'temperature', 'replicate', parameter)]
  # remove potential percentage signs and make value numeric
  QC[,parameter] <- as.numeric(sub('%', '', QC[,parameter]))
  QC$acnTemp <- paste(QC$accession, QC$temperature, sep = '_')
  QC$acnTemp <- factor(QC$acnTemp, ordered = TRUE)
  QC$tempCol <- as.character(tempCol$col[match(QC$temperature, tempCol$temp)])
  QC$repSymbol <- repSymb$symb[match(QC$replicate, repSymb$rep)]
  # set graphical parameters to allow for legdn outside the plot
  opar <- par
  par(xpd = T, mar = par()$mar + c(0,0,0,4))
  plot(as.numeric(QC$acnTemp), QC[,parameter], type = 'p', pch = QC$repSymbol, col = QC$tempCol, ylab = parameter, xaxt = 'n', xlab = '', main = parameter)
  tempLabels <- unlist(lapply(sort(unique(QC$acnTemp)), FUN = function(x){strsplit(as.character(x), '_')[[1]][2]}))
  tempTicks <- unique(as.numeric(sort(QC$acnTemp)))
  acnLabels <- unique(unlist(lapply(sort(unique(QC$acnTemp)), FUN = function(x){strsplit(as.character(x), '_')[[1]][1]})))
  acnTicks <- tempTicks[seq(1, length(tempTicks), 2)] + 0.5

  axis(1, at = tempTicks, labels = tempLabels)
  axis(1, at = acnTicks, labels = acnLabels, line = 1, tick = F, col = NA)
  
  legendNames <- c(as.character(temperatures), as.character(replicates))
  legendColours <- c(as.character(tempCol$col), rep('black', length(replicates)))
  legendLwd <- c(rep(4, length(temperatures)), rep(NA, length(replicates)))
  legendSymbols <- c(rep(NA, length(temperatures)), repSymb$symb)

  legend(max(tempTicks) +0.5, max(QC[,parameter]), legend = legendNames, col = legendColours, pch = legendSymbols, lwd = legendLwd, bty = 'n')
  par <- opar
}

#PCA with colours according temperature and shape according accession
PCAplot_temp_acn <- function(dds)
{
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("temperature", "accession"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=accession)) +
    geom_point(size=3) +
    scale_color_manual(values=c('#0047BB', '#8BD52D')) +
    scale_shape_manual(values=c(15, 19, 17, 18, 4, 8, 9, 11)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
}

#PCA with colours according temperature and shape according replicate
PCAplot_temp_rep <- function(dds, title = '')
{
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("temperature", "replicate"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=replicate)) +
    geom_point(size=3) +
    scale_color_manual(values=c('#0047BB', '#8BD52D')) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(title) +
    coord_fixed()
}



# function to create dataframe with differentially expressed genes selected on log2FoldChange and adjusted p-value
# ready to write to file
getDF_select_FC_padj <- function(res.deseq, FC = 2, padj = 0.05)
{
  res.deseq <- na.omit(res.deseq)
  res.deseq.FC_padj <- res.deseq[abs(res.deseq$log2FoldChange) > FC & res.deseq$padj < padj, ]
  summary(res.deseq.FC_padj)
  res.deseq.FC_padj <- as.data.frame(res.deseq.FC_padj)
  res.deseq.FC_padj$gene <- row.names(res.deseq.FC_padj)
  res.deseq.FC_padj <- res.deseq.FC_padj[, c(7,1,2,3,4,5,6)]
  return(res.deseq.FC_padj)
}


# volcanoplot
volcanoplot_geneSelection <- function(genes, background, main = 'volcano plot', FC = 2, padj = 0.05)
{
  if (!is.data.frame(background))
  {
    background <- na.omit(as.data.frame(background))
  }
  genes.data <- background[genes, ]
  plot(background$log2FoldChange, -log10(background$padj), pch=20, col = 'grey', main = main, cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value))
  points(genes.data$log2FoldChange, -log10(genes.data$padj), pch=20, col = '#0047BB', cex=1.0)
  
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-FC, col="black", lty=4, lwd=2.0)
  abline(v=FC, col="black", lty=4, lwd=2.0)
  abline(h=-log10(padj), col = 'black', lty = 4, lwd = 2.0)
}

# plot normalised counts of specific genes
plotNormCounts <- function(dds, genes, samples, legend = F)
{
  genes <- as.character(genes)
  dds <- estimateSizeFactors(dds)
  dds.normCounts <- counts(dds, normalized=TRUE)
  
  normCounts <- as.data.frame(dds.normCounts)
  # select genes
  normCounts <- normCounts[genes, ]
  
  # log transform
  normCounts  <- log(normCounts)
  minCount <- 0
  maxCount <- max(normCounts, na.rm = T)
  
  # make legend
  geneCol <- 'black'
  if (legend == TRUE)
  {
    geneCol <- rainbow(n= length(genes))
    
  }
  
  # prepare dataframe
  normCounts <-  reshape(normCounts, direction = "long", varying = names(normCounts), v.names = "normCount", 
                         ids = rownames(normCounts), timevar = "sample", times = names(normCounts))
  normCounts$temperature <- samples$temperature[match(normCounts$sample, samples$sample)]
  normCounts$temperature <- factor(normCounts$temperature, levels = c('6C', '16C'), ordered = T)
  normCounts$accession <- samples$accession[match(normCounts$sample, samples$sample)]
  avgNormCounts <- aggregate(normCount ~ temperature * id, data = normCounts, FUN = 'mean')
  avgNormCounts$temperature <- factor(avgNormCounts$temperature, levels = c('6C', '16C'), ordered = T)
  
  plot(NA, NA, xlim = c(0.5, 2.5), ylim = c(minCount, maxCount), ylab = 'log(normalised counts)', xaxt = 'n', xlab = '')
  axis(side = 1, at = c(1,2), labels = c('6C', '16C'))
  for (g in genes)
  {
    avgNormCounts.sub <- avgNormCounts[avgNormCounts$id == g, ]
    
    lines(avgNormCounts.sub$temperature, avgNormCounts.sub$normCount, col = geneCol[match(g, genes)])
  }
  legend('topright', legend = genes, col = geneCol, lwd = rep(2, length(genes)), bty = 'n')
}

# check for two genes if the log2foldchange have opposite directions
checkForOppositeDiffExpr <- function(geneX, geneY, res, pvalThrs)
{
  trueOpposites <- FALSE
  lfcX <- res[geneX, 'log2FoldChange']
  padjX <- res[geneX, 'padj']
  lfcY <- res[geneY, 'log2FoldChange']
  padjY <- res[geneY, 'padj']
  if (padjX > pvalThrs | is.na(padjX)){lfcX <- 0}
  if (padjY > pvalThrs | is.na(padjY)){lfcY <- 0}
  if ((lfcX > 0 & lfcY < 0) | (lfcX < 0 & lfcY > 0))
  {
    trueOpposites <- TRUE
  }
  return(trueOpposites)
}

# check for two genes if the log2foldchange cross each other
checkForCrossedDiffExpr <- function(geneX, geneY, cts)
{
  cross <- FALSE
  d6C_XY <- cts$vst6C[cts$gene == geneX] - cts$vst6C[cts$gene == geneY]
  d16C_XY <- cts$vst16C[cts$gene == geneX] - cts$vst16C[cts$gene == geneY]
  
  if ((d6C_XY > 0 & d16C_XY < 0) | (d6C_XY < 0 & d16C_XY > 0))
  {
    cross <- TRUE
  }
  return(cross)
}

# do a sleection on deseq results for specific fold change and adjusted p-value cut-off
getResults <- function(res, FC, padj)
{
  res <- na.omit(res)
  res_FC.padj <- res[res$padj < padj & abs(res$log2FoldChange) >FC, ]
  return(res_FC.padj)
}


