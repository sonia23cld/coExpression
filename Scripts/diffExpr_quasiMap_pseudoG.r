# differential expresison analysis from quasiMapped (salmon) data based pseudoGenomes and Araport11 annotation


# LIBRARIES #
library(tximport)
library(DESeq2)
library(tidyverse)
library(vsn)
library(ggplot2)
source('/Volumes/nordborg/user/sonia.celestini/coExpression/Scripts/diffExpr_functions.r')

# DATA #
#salmon_qc
#salmon_qc <- read.csv('Results/SalmonQuant/salmon_quasiMap_pseudoG_summary.csv')
#salmon_qc <- salmon_qc[!salmon_qc$temperature %in% c('16C?', '6C?'), ]

# QC #
# quality check the quasi mapping - mapping rates
#pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/mappingRate_salmon_quasiMap_pseudoG.pdf')
#QCplot(salmon_qc, 'mappingRate')
#dev.off()

#import sample data
samples <- read.table('/Volumes/nordborg/user/sonia.celestini/WGCNA_8acn_new/samples_updated.txt', header = T, sep = ' ', comment.char = "")
# filter only sequenced samples
samples <- filter(samples, Selected == 'yes') %>%
  mutate(accession = as.factor(accession),
         temperature = as.factor(temperature),
         replicate = as.factor(replicate),
         experiment = as.factor(paste(temperature, replicate, sep = '_')),
         ID = paste(tray, tray_coordinate, temperature, replicate, sep = '_'),
         sampleName = paste(accession, temperature, replicate, sep = '_'))

# file translating transcipts to genes
Araport11 <- read.table('/Volumes/nordborg/user/sonia.celestini/WGCNA_8acn_new/Araport11_GFF3_genes_transposons.201606.gtf', sep ='\t') 

# prepare Araport11 data
colnames(Araport11) <- c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
Araport11$attributes <- as.character(Araport11$attributes)
library(stringr)
Araport11[,9:10]<- str_split_fixed(Araport11$attributes, ';', 3)[,1:2]
Araport11$V10<- sub(".*id ", "", Araport11$V10)

# create transcript to gene annotation dataframe
tx2gene <- data.frame('transcriptID' = sub(".*id ", "", Araport11$attributes), 'geneID' = Araport11$V10)
# data
expression.files <- file.path('/Volumes/nordborg/pub/forPieter/WGCNA/WGCNA_8acn/SalmonQuant/', samples$basename,'_quasiMap_pseudoG_Trimmed/quant.sf', fsep = '')
names(expression.files) <- samples$sample

# import salmon quantification data
samples.acn1 <- as.character(samples$sample[samples$accession %in% unique(samples$accession)[1]])
txi <- tximport(expression.files[samples.acn1], tx2gene = tx2gene, type = 'salmon', dropInfReps = T)
for (acn in unique(samples$accession)[-1]) {
  samples.acn <- as.character(samples$sample[samples$accession %in% acn])
  txi.acn <- tximport(expression.files[samples.acn], tx2gene = tx2gene, type = 'salmon', dropInfReps = T)
  # use overlapping genes
  idx <- intersect(rownames(txi$counts), rownames(txi.acn$counts))
  txi$abundance <- cbind(txi$abundance[idx, ], txi.acn$abundance[idx, ])
  txi$counts <- cbind(txi$counts[idx, ], txi.acn$counts[idx, ])
  txi$length <- cbind(txi$length[idx, ], txi.acn$length[idx, ])
}
#match Colnames with samples order
txi$abundance <- txi$abundance[ ,as.vector(samples$sample)]
txi$counts <- txi$counts[ ,as.vector(samples$sample)]
txi$length <- txi$length[ ,as.vector(samples$sample)]

#DESeq
dds.full <- DESeqDataSetFromTximport(txi, colData = samples, design = ~  replicate + accession + temperature + replicate:temperature + accession:temperature)
dds.temp <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ replicate + accession + temperature + replicate:temperature)

# PCA quantification data
#pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/PCA/PCA_temp_acn_quasiMap_pseudoG.pdf')
PCAplot_temp_acn(dds.full)
dev.off()

# PCA plots per accession
library(cowplot)
p <- list()
for (acn in unique(samples$accession)) {
  samples.acn <- samples[samples$accession == acn, ]
  txi.acn <- tximport(expression.files[as.vector(samples.acn$sample)], tx2gene = tx2gene, type = 'salmon', dropInfReps = T)
  dds.a.full <- DESeqDataSetFromTximport(txi.acn, colData = samples.acn, design = ~  replicate + temperature)
  #pdf(paste('Results/Transcriptome/DiffExpr_allAcns/Plots/PCA/PCA_quasiMap_pseudoG_', a, '.pdf', sep = ''))
  vsd <- vst(dds.a.full, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("temperature", "replicate"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p[[acn]]<- ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=replicate)) +
    geom_point(size=3) +
    scale_color_manual(values=c('#0047BB', '#8BD52D')) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(acn) +
    theme(aspect.ratio=1, legend.position = "none")
  #dev.off()
}

legend <- cowplot::get_legend(p[[1]] + theme(legend.position = "bottom"))
p_grid <- cowplot::plot_grid(plotlist = p, ncol = 4)
cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(2, 0.2))

save(dds.full, file = '/Volumes/nordborg/user/sonia.celestini/WGCNA_8acn_new/dds_full.Rdata')
#save(dds.temp, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/dds.temp.quasiMap.pseudoG.Rdata')

# DE analysis over all accessions
# set reference condition
dds.full$temperature <- relevel(dds.full$temperature, ref = "16C")
dds.temp$temperature <- relevel(dds.temp$temperature, ref = "16C")
dds.full$accession <- relevel(dds.full$accession, ref = '6017')
dds.temp$accession <- relevel(dds.temp$accession, ref = '6017')


# differential expression analysis
dds.full.deseq <- DESeq(dds.full)
dds.temp.deseq <- DESeq(dds.temp)

# Main effect
# 6vs 16 over all accessions
res.6vs16 <- results(dds.temp.deseq, contrast = c("temperature", "6C", "16C"))
res.6vs16.FC2_p5 <- getDF_select_FC_padj(res.6vs16, FC = 2, padj = 0.05)
save(res.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.6vs16.quasiMap.pseudoG.Rdata')
write.table(res.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.6vs16.FC2_p5$gene, res.6vs16, main = '6vs16', FC = 2, padj = 0.05)
dev.off()


# 6vs16 for each separate accession
# the temperature effect for 6017 (the main effect)
res.6017.6vs16 <- results(dds.full.deseq, contrast=c("temperature","6C","16C"))
res.6017.6vs16.FC2_p5 <- getDF_select_FC_padj(res.6017.6vs16, FC = 2, padj = 0.05)
#save(res.6017.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.6017.6vs16.quasiMap.pseudoG.Rdata')
#write.table(res.6017.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_6017_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)

#temperature effect for the rest of the accessions
res.6Vs16_allacn<- list()
res.6Vs16_allacn$`6017`$result<- res.6017.6vs16
res.6Vs16_allacn$`6017`$FC2_p5<- res.6017.6vs16.FC2_p5

for (acn in unique(samples$accession)[-1]) {
  res.6Vs16_allacn[[acn]]$result<- results(dds.full.deseq, contrast=list( c("temperature_6C_vs_16C", paste("accession", acn, ".temperature6C", sep = ''))))
  res.6Vs16_allacn[[acn]]$FC2_p5<- getDF_select_FC_padj(res.6Vs16_allacn[[acn]]$result, FC = 2, padj = 0.05)
}

#plot volcano plots of the results for all the accessions together
par(mfrow=c(2, 4), pty='s')
for (acn in unique(samples$accession)) {
  volcanoplot_geneSelection(res.6Vs16_allacn[[acn]]$FC2_p5$gene,  res.6Vs16_allacn[[acn]]$result, main = paste('6vs16 in ', acn, sep=''), FC = 2, padj = 0.05)
}

#investigate differential expression (AllDE, UpDE, DownDE)
sumtable<-data.frame('accession'= unique(samples$accession), 'AllDE'=NA, 'UpDE'=NA, 'DownDE'=NA)
#loop creation to fill the table 
for (acn in unique(samples$accession)){
  sumtable[sumtable$accession==acn,'AllDE']<- length(res.6Vs16_allacn[[acn]]$FC2_p5$gene)
  acn.data.up<-res.6Vs16_allacn[[acn]]$FC2_p5[res.6Vs16_allacn[[acn]]$FC2_p5$log2FoldChange>0,]
  sumtable[sumtable$accession==acn,'UpDE']<-length(acn.data.up$gene)
  acn.data.down<-res.6Vs16_allacn[[acn]]$FC2_p5[res.6Vs16_allacn[[acn]]$FC2_p5$log2FoldChange<0,]
  sumtable[sumtable$accession==acn,'DownDE']<-length(acn.data.down$gene)
}
write.table(sumtable,file='/Volumes/nordborg/user/sonia.celestini/WGCNA_8acn_new/Gene_DiffExpr.txt',quote = FALSE, row.names = FALSE, col.names = T, sep = '\t')


#=======================================
#I didn't use this part
# LRT to find genes with accession specific responses
dds.acnSpec <- estimateSizeFactors(dds.full)
dds.acnSpec <- estimateDispersions(dds.acnSpec)
dds.acnSpec <- nbinomLRT(dds.acnSpec, full = ~ replicate + accession + temperature + accession:temperature, reduced = ~ replicate + accession + temperature)
res.acnSpec <- results(dds.acnSpec)
res.acnSpec.FC0_p10 <- getDF_select_FC_padj(res.acnSpec, FC = 0, padj = 0.1)
save(res.acnSpec, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.acnSpec.quasiMap.pseudoG.Rdata')
write.table(res.acnSpec.FC0_p10, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_acnSpec_FC0_p10_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
res.acnSpec.FC2_p5 <- getDF_select_FC_padj(res.acnSpec, FC = 2, padj = 0.05)
write.table(res.acnSpec.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_acnSpec_FC2_p5_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16 accession specific
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_acnSpec_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.acnSpec.FC0_p10$gene, res.acnSpec, main = ' accession specific 6vs16', FC = 0, padj = 0.1)
dev.off()

# the interaction term for condition effect in 6909 vs 6017.
# this tests if the temperature effect is different in 6909 compared to 6017
res.6vs16C.6909vs6017 <- results(dds.deseq, name="accession6909.temperature6C")
# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))
# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.





