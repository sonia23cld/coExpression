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
  output <- data.frame(name = I(array("",dim=c(1,n_res))),
                       log10_pval = I(array(NA,dim=c(1,n_res))),
                       group_index = I(array(NA,dim=c(1,n_res))),
                       groups_size = I(array(NA,dim=c(1,n_res))),
                       sample_size = I(array(NA,dim=c(1,n_res))),
                       intersect_size = I(array(NA,dim=c(1,n_res))),
                       set = I(array(NA,dim=c(1,n_res))))
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
Enrichment<- AraEnrice_v1(genes, OG, n_res = 50)

#repeat only for GxT genes in module 31
Gene_info<- read.table('/Volumes/nordborg/user/sonia.celestini/WGCNA_newData/Gene_category.txt', header = T)
GxT<- Gene_info$genes[Gene_info$category %in% 'GxT']
gxt_31<- Genes_31[Genes_31 %in% GxT]

genes_gxt<- c()
for(i in 1:length(gxt_31)) {    
  genes_gxt <- c(genes_gxt, which(gene_infoV2$Name == gxt_31[i]))
}

Enrichment_gxt<- AraEnrice_v1(genes_gxt, OG, n_res = 50)

