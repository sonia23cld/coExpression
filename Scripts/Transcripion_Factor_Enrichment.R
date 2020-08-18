#transcription factors enrichment for network modules
##used script from Yoav

library(Matrix)

#define the function

AraEnrice_v1 <- function(genes, OG, bg_genes = NULL, n_res = 25) {
  genes2vec <- function(genes, n) {
    vec <- array(0, dim = c(1,n))
    vec[genes] <- 1  
    return(vec)
  }
  n_genes <- dim(OG$groups)[1];
  #############################################################################
  # If no background genes where specified, all genes will be considerd
  if(is.null(bg_genes)) {bg_genes <- 1:n_genes}
  bg_vec <- genes2vec(bg_genes, n_genes)
  #############################################################################
  
  #############################################################################
  # We want to caluclate a hypergeometric test on our genes vs. all the groups
  genes_vec <- genes2vec(intersect(genes, bg_genes), n_genes)
  # Size of the spaces - information for every group
  space_size <- bg_vec %*% OG$groups_bg
  space_size <- space_size[OG$g_info$background_index] 
  # size of sample genes
  sample_size <- genes_vec %*% OG$groups_bg
  sample_size <- sample_size[OG$g_info$background_index] 
  # size of groups
  groups_size <-  as.array(bg_vec %*% OG$groups)
  # size of intersect (interesting part)
  intersect_size <-  as.array(genes_vec %*% OG$groups)
  
  # Calculate the hypergeometric test
  res <- phyper(q=intersect_size-1,
                m = groups_size,
                n = space_size-groups_size,
                k = sample_size,
                log.p = TRUE, lower.tail = FALSE) / log(10)
  #############################################################################
  # Plot the results!
  sort_res <- sort(res, index.return = TRUE, decreasing = FALSE)
  output <- data.frame(matrix(NA, nrow = n_res, ncol = 8)) 
  colnames(output)<- c('name', 'log10_pval', 'padj', 'group_index', 'groups_size', 'sample_size', 'intersect_size', 'set')    
  writeLines("using log10") 
  for(i in 1:n_res) {
    g_i <- sort_res$ix[i]
    output$group_index[i] <- g_i
    output$name[i] <- OG$g_info$name[g_i]
    output$log10_pval[i] <- round(sort_res$x[i],digits = 1)
    output$groups_size[i] <- groups_size[g_i]
    output$intersect_size[i] <- intersect_size[g_i]
    output$sample_size[i] <- sample_size[g_i]
    output$set[i] <- OG$g_info$set[g_i]
    
    rm(list = c('g_i'))
    
    writeLines(paste(
      i, "\t[", output$log10_pval[i] , "] [", 
      output$groups_size[i]," | ", output$sample_size[i], " | ", output$intersect_size[i], "]\t",
      output$set[i],':',output$name[i],'\r\n',
      sep=""))
  }
  output$padj<- p.adjust(10^output$log10_pval, method = 'BH', n=349)
  invisible(output)
}

#=================
#start with my data
load("/Volumes/nordborg/pub/1001T/organize_groups_20191002.Rdata") 
load("/Volumes/nordborg/pub/1001T/gene_infoV2.Rdata")
#load genes you want to analyze, in my case genes from module 31
load('/Volumes/nordborg/user/sonia.celestini/WGCNA_newData/Network.RData')
Genes_31<- colnames(all.expr$data)[all.expr$mergedLabels == 31]

#to find enrichment against binding of TFs you need data from DAP-seq. You can see it's present in OG typing:
unique(OG$g_info$set)

#create vector of genes indexes
genes<- c()

for(i in 1:length(Genes_31)) {    
  genes <- c(genes, which(gene_infoV2$Name == Genes_31[i]))
}

#I want results only for DAP-seq, so I have to change OG
ind <- which(OG$g_info$set == "DAP-seq applied on DNA with modifications")
length(ind)
OG$g_info <- OG$g_info[ind,]
OG$groups <- OG$groups[,ind]

#TF Enrichment
Enrichment<- AraEnrice_v1(genes, OG, n_res = 349)
Enrichment<- Enrichment[Enrichment$padj < 0.05,]
write.table(Enrichment, file = '/Volumes/nordborg/pub/forPieter/WGCNA/WGCNA_8acn/TF_enrichment_module31.txt', sep = '\t', quote = F, row.names = F)

#repeat only for GxT genes in module 31
Gene_info<- read.table('/Volumes/nordborg/user/sonia.celestini/WGCNA_newData/Gene_category.txt', header = T)
GxT<- Gene_info$genes[Gene_info$category %in% 'GxT']
gxt_31<- Genes_31[Genes_31 %in% GxT]

genes_gxt<- c()
for(i in 1:length(gxt_31)) {    
  genes_gxt <- c(genes_gxt, which(gene_infoV2$Name == gxt_31[i]))
}

Enrichment_gxt<- AraEnrice_v1(genes_gxt, OG, n_res = 349)
Enrichment_gxt<- Enrichment_gxt[Enrichment_gxt$padj < 0.05,]
write.table(Enrichment_gxt, file = '/Volumes/nordborg/pub/forPieter/WGCNA/WGCNA_8acn/TF_enrichment_module31_onlyGxTgenes.txt', sep = '\t', quote = F, row.names = F)

###########################################
#check which genes have a specific TF motif
#check for CBF family
check_Tf<-OG$groups[genes, Enrichment$group_index[c(1,7, 11, 12)]] #here you mean whether your genes are in group 23 ecc...? I.e. the syntax is OG$groups[genei, groupj]?
#If they are present you will have a 1 

CBF_genes<-which(apply(check_Tf, 1, function(x) {any(x %in% 1)}))
CBF_genes<- genes[CBF_genes]
CBF_genes<-  gene_infoV2$Name[CBF_genes]

#check for DDF family
check_Tf<-OG$groups[genes, Enrichment$group_index[c(3, 6)]] 
DDF_genes<-which(apply(check_Tf, 1, function(x) {any(x %in% 1)}))
DDF_genes<- genes[DDF_genes]
DDF_genes<-  gene_infoV2$Name[DDF_genes]

###########################################

#you will have to do SNPs analysis on significant TF genes you spotted. You will get the gene name coding for the TF in the gff file.
TF_families<- unique(c(Enrichment$name, Enrichment_gxt$name))

#impot and clean gff file
Araport11 <- read.gff('/Volumes/nordborg/user/sonia.celestini/Araport11_GFF3_genes_transposons.201606.gff') 
Araport11[,9:10]<- str_split_fixed(Araport11$attributes, ';', 2)
Araport11<- Araport11[Araport11$type %in% 'gene',]
Araport11$attributes<- str_split_fixed(Araport11$attributes, 'ID=', 2)[,2]
colnames(Araport11)[9]<- 'ID'
colnames(Araport11)[10]<- 'attributes'


#clean TF names
TF_families<- str_split_fixed(TF_families, '_', 2)[,1]

#get the genes name
TF_families<- gsub("t", "T", TF_families)
TF_families<- gsub("g", "G", TF_families)
TF_genes<- list()
for (i in TF_families) {
  TF_genes[[i]]<- Araport11$ID[grep(i, Araport11$attributes)]
}

save(TF_genes, file = '/Volumes/nordborg/user/sonia.celestini/WGCNA_newData/TFgenes_module31.RData')





