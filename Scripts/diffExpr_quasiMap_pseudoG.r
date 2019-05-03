# differential expresison analysis from quasiMapped (salmon) data based pseudoGenomes and Araport11 annotation

setwd('/Volumes/nordborg/user/pieter.clauw/Documents/Experiments/UltimateQandD/')

# LIBRARIES #
library(tximport)
library(DESeq2)
library(vsn)
library(ggplot2)
source('Scripts/Transcriptome/diffExpr_functions.r')

# DATA #
#salmon_qc
salmon_qc <- read.csv('Results/Transcriptome/SalmonQuant/Araport11/salmon_quasiMap_pseudoG_summary.csv')
salmon_qc <- salmon_qc[!salmon_qc$temperature %in% c('16C?', '6C?'), ]

#update samples accordingly, original file in projects folder
samples <- read.table('Data/Transcriptome/samples.txt', header = T, sep = '\t', comment.char = "")
samples$basename <- as.character(samples$basename)
samples$accession <- as.factor(samples$accession)
samples$acnTemp <- paste(samples$accession, samples$temperature, sep = '_')
row.names(samples) <- samples$sample

# only samples with confirmed genotypes
samples <- samples[samples$genoCheck == 'ok', ]
# remove samples 63076 and 63077 -> unlcear which is what temperature, likely to be swapped
samples <- samples[!samples$sample %in% c(63076, 63077), ]
samples$temperature <- factor(as.character(samples$temperature))
table(samples[,c('accession', 'temperature', 'replicate')])

files <- file.path('Results/Transcriptome/SalmonQuant/Araport11/', samples$basename, '' ,'_quasiMap_pseudoG/quant.sf', fsep = '')
names(files) <- samples$sample

# file translating transcipts to genes
Araport11 <- read.table('/Volumes/nordborg/user/pieter.clauw/Documents/Source/Araport11/Araport11_GFF3_genes_transposons.201606.ChrM_ChrC_FullName.gtf')
colnames(Araport11) <- c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
Araport11$attributes <- as.character(Araport11$attributes)

tx2gene <- data.frame('transcriptID' = Araport11[,10], 'geneID' = Araport11[,13])

# import salmon quantification data
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, dropInfReps = T)
dds.full <- DESeqDataSetFromTximport(txi, colData = samples, design = ~  replicate + accession + temperature + accession:temperature)
dds.temp <- DESeqDataSetFromTximport(txi, colData = samples, design = ~  replicate + accession + temperature)

save(dds.full, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/dds.full.quasiMap.pseudoG.Rdata')
save(dds.temp, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/dds.temp.quasiMap.pseudoG.Rdata')

# QC #
# quality check the quasi mapping - mapping rates
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/mappingRate_salmon_quasiMap_pseudoG.pdf')
QCplot(salmon_qc, 'mappingRate')
dev.off()

# PCA quantification data
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/PCA/PCA_temp_acn_quasiMap_pseudoG.pdf')
PCAplot_temp_acn(dds.full)
dev.off()

# PCA plots per accession
for (a in unique(samples$accession))
{
  samples.a <- samples[samples$accession == a, ]
  files.a <- file.path('Results/Transcriptome/SalmonQuant/Araport11/', samples.a$basename, '' ,'_quasiMap_pseudoG/quant.sf', fsep = '')
  names(files.a) <- samples.a$sample
  txi.a <- tximport(files.a, type = 'salmon', tx2gene = tx2gene, dropInfReps = T)
  dds.a.full <- DESeqDataSetFromTximport(txi.a, colData = samples.a, design = ~  replicate + temperature)
  pdf(paste('Results/Transcriptome/DiffExpr_allAcns/Plots/PCA/PCA_quasiMap_pseudoG_', a, '.pdf', sep = ''))
  print(PCAplot_temp_rep(dds.a.full, title = a))
  dev.off()
}

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
save(res.6017.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.6017.6vs16.quasiMap.pseudoG.Rdata')
write.table(res.6017.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_6017_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16 in 6017
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_6017_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.6017.6vs16.FC2_p5$gene, res.6017.6vs16, main = '6vs16 in 6017', FC = 2, padj = 0.05)
dev.off()

# the temperature effect for 6909.
# this is the main effect *plus* the interaction term
# (the extra condition effect in 6909 compared to 6017).
res.6909.6vs16 <- results(dds.full.deseq, contrast=list( c("temperature_6C_vs_16C","accession6909.temperature6C") ))
res.6909.6vs16.FC2_p5 <- getDF_select_FC_padj(res.6909.6vs16, FC = 2, padj = 0.05)
save(res.6909.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.6909.6vs16.quasiMap.pseudoG.Rdata')
write.table(res.6909.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_6909_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16 in 6909
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_6909_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.6909.6vs16.FC2_p5$gene, res.6909.6vs16, main = '6vs16 in 6909', FC = 2, padj = 0.05)
dev.off()

# the temperature effect for 9559
# this is the main effect *plus* the interaction term
# (the extra condition effect in 9559 compared to 6017).
res.9559.6vs16 <- results(dds.full.deseq, contrast=list( c("temperature_6C_vs_16C","accession9559.temperature6C") ))
res.9559.6vs16.FC2_p5 <- getDF_select_FC_padj(res.9559.6vs16, FC = 2, padj = 0.05)
save(res.9559.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.9559.6vs16.quasiMap.pseudoG.Rdata')
write.table(res.9559.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_9559_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16 in 9559
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_9559_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.9559.6vs16.FC2_p5$gene, res.9559.6vs16, main = '6vs16 in 9559', FC = 2, padj = 0.05)
dev.off()

# the temperature effect for 9728
# this is the main effect *plus* the interaction term
# (the extra condition effect in 9728 compared to 6017).
res.9728.6vs16 <- results(dds.full.deseq, contrast=list( c("temperature_6C_vs_16C","accession9728.temperature6C") ))
res.9728.6vs16.FC2_p5 <- getDF_select_FC_padj(res.9728.6vs16, FC = 2, padj = 0.05)
save(res.9728.6vs16, file = 'Results/Transcriptome/DiffExpr_allAcns/DESeq2_results/res.9728.6vs16.quasiMap.pseudoG.Rdata')
write.table(res.9728.6vs16.FC2_p5, file = 'Results/Transcriptome/DiffExpr_allAcns/6vs16_9728_quasiMap_pseudoG.txt', sep = '\t', row.names = F, quote = F)
# voclanoPlot 6vs16 in 9728
pdf('Results/Transcriptome/DiffExpr_allAcns/Plots/VolcanoPlots/volcano_6vs16_9728_quasiMap_pseudoG.pdf')
volcanoplot_geneSelection(res.9728.6vs16.FC2_p5$gene, res.9728.6vs16, main = '6vs16 in 9728', FC = 2, padj = 0.05)
dev.off()


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





