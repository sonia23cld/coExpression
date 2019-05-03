#Analyse overlap
#create graphs of modules expression behavior and find Go terms in common

#find Go terms in common 
#before, you have to do the enrichment of the modules you want
#use annotation and geneuniverse saved
#ClearTable2 represent a table with the module names you want 
library(topGO)
for (i in 1:nrow(ClearTable2)) {
    #extract genes (it doesn' matter which table you use, they are the same genes)
    modulename<- ClearTable2$NameModule6[[i]]
    mygenes<-colnames(multiExpr[[2]]$data)[moduleColors6C == modulename]
    
    # remove any candidate genes without GO annotation
    keep = mygenes %in% go_ids[,2]
    keep =which(keep==TRUE)
    mygenes=mygenes[keep]
    
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
    
    #get list of significant GO before multiple testing correction
    results.table.p= all_res_final[which(all_res_final$weightFisher<=0.001),]
    
    #get list of significant GO after multiple testing correction
    results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
    
    #create a specific folder for the results
    dir.create(paste('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/Go enrichment6C/', modulename, sep = ''))
    
    #save first top 50 ontolgies sorted by adjusted pvalues
    write.table(all_res_final[1:50,], paste('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/Go enrichment6C/', modulename, '/summary_topGO_analysis_', modulename, '.csv', sep=''), sep=",",quote=FALSE,row.names=FALSE)
    
    # PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
    
    pdf(file= paste('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/Go enrichment6C/', modulename, '/topGOPlot_fullnames_', modulename, '.pdf', sep = ''), height=12, width=12, paper='special', pointsize=18)
    showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "all", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
    dev.off()
    #get genes in significant GO terms
    myterms =results.table.p$GO.ID 
    if (!length(myterms) == 0) {
      mygenes = genesInTerm(GOdata, myterms)
      
      var=c()
      for (y in 1:length(myterms)) {
        myterm=myterms[y]
        mygenesforterm= mygenes[myterm][[1]]
        mygenesforterm=paste(mygenesforterm, collapse=',')  
        var[y]=paste("GOTerm",myterm,"genes-",mygenesforterm)
      }
      
      write.table(var, paste('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/Go enrichment6C/', modulename, '/genetoGOmapping_', modulename, '.txt', sep = ''), sep="\t",quote=F)
    }
} 

#Goterms in common
library(officer)
setwd('/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/')

#trial with doc (not finished)
doc<-read_docx()
doc<- doc %>% 
  body_add_par(value = 'GO terms overlapping between 16C and 6C', style = 'centered')

#trial with table
TableGO<- data.frame('Name module 16C'= as.character(), 'Name module 6C'= as.character(), 'GO term'= as.character(), 'Function'= as.character())
for (j in 1:nrow(ClearTable2)) {
  namemodule6<- as.character(ClearTable2$NameModule6[[j]])
  namemodule16<- as.character(ClearTable2$NameModule16[[j]])
  table6<- read.table(paste('Go enrichment6C/', namemodule6, '/summary_topGO_analysis_', namemodule6, '.csv', sep=''),  header = F, sep = ',', fill = TRUE, as.is = T)
  table16<- read.table(paste('Go enrichment16C/', namemodule16, '/summary_topGO_analysis_', namemodule16, '.csv', sep=''),  header = F, sep = ',', fill = TRUE, as.is = T, quote = '')
  colnames(table6)<- table6[1,]
  colnames(table16)<- table16[1,]
  table6<-table6[-1,]
  table16<-table16[-1,]
  table6<- table6[,-8]
  table16<- table16[,-8]
  table6<-table6[table6$GO.ID!=1,]
  table16<-table16[table16$GO.ID!=1,]
  for (x in 1:length(table6$GO.ID)) {
    GOname6<-table6$GO.ID[[x]]
    for (y in 1:length(table16$GO.ID)) {
      GOname16<- table16$GO.ID[[y]]
      if (GOname6 == GOname16) {
        lineiwant<- data.frame('Name module 16C'= namemodule16, 'Name module 6C'= namemodule6, 'GO term'= GOname16, 'Function'= table16$Term[[y]])
        TableGO<-rbind(TableGO, lineiwant)
      next}
    }
  }
}

write.table(TableGO, file= '/Volumes/nordborg/pub/forPieter/WGCNA/Results/Results_Mendel/New_parameters_soft14/Table GOterms in common 16C vs 6C.txt', quote = FALSE, sep = '\t', col.names = TRUE)



#behaviour of the expression for the Goterms in common 

#obtain genes in the Goterm (use geneuniverse and annotation saved)
library(topGO)
# make named factor showing which genes are of interest
geneList=factor(as.integer(genesuniverse %in% colnames(multiExpr[[1]]$data)))
names(geneList)= genesuniverse
#make topGO data object
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
#obtain genes
mygenes = genesInTerm(GOdata, Goterms)

#create graph of behaviour
#create a table with expression of genes of interest for a certain Go term
Tableplot<- rbind(multiExpr[[1]]$data[, colnames(multiExpr[[1]]$data) %in% mygenes], multiExpr[[2]]$data[, colnames(multiExpr[[2]]$data) %in% mygenes])

#define order (I took 'sample' from WGCNA script)
df<- samples[,1:3]
df6<-df[df$temperature %in% '6C',]
df6<-df6[order(df6$accession),]
df16<-df[df$temperature %in% '16C', ]
df16<-df16[order(df16$accession),]
orderiwant<-c(df16$sample, df6$sample)
orderiwant<-as.factor(orderiwant)

#order the table
Tableplot<-Tableplot[match(orderiwant, row.names(Tableplot)),]

#get the range for x and y axis
yrange<- range(Tableplot)
xrange<- range(1:nrow(Tableplot))

#set up the plot
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plot(NA, NA, xlim=xrange, ylim = yrange, xaxt = 'n', type='n', xlab='Name sample', ylab='Expression level') #main=paste('Expression level genes in module consensus brown', Goterms[[1]], sep = ''))
axis(1, at=1:22, labels = rownames(Tableplot), cex.axis=0.5)
colors <- rainbow(ncol(Tableplot))
legend('topright', inset = c(-0.2,-0.1), colnames(Tableplot), fill =colors, cex=0.5)
#add lines
for (i in 1:ncol(Tableplot)) {
  lines(Tableplot[i], type='l', lwd=1.5, col=colors[[i]])
}
dev.off() 

#genes for module consensus 
mygenes<-colnames(multiExpr[[1]]$data)[moduleColors== 'sienna2']
Tableplot<- rbind(multiExpr[[1]]$data[, colnames(multiExpr[[1]]$data) %in% mygenes], multiExpr[[2]]$data[, colnames(multiExpr[[2]]$data) %in% mygenes])
Tableplot<-Tableplot[match(orderiwant, row.names(Tableplot)),]
yrange<- range(Tableplot)
xrange<- range(1:nrow(Tableplot))

#set up the plot
groups<-c('16C', '6C')
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
plot(NA, NA, xlim=xrange, ylim = yrange, xaxt = 'n', type='n', xlab='Name sample', ylab='Expression level') #main=paste('Expression level genes in module consensus brown', Goterms[[1]], sep = ''))
axis(1, at=1:22, labels = rownames(Tableplot), cex.axis=0.5)
colors <- rainbow(ncol(Tableplot))
abline(v=12.5, lty=2)
#legend('topright', inset = c(-0.2,-0.1), colnames(Tableplot), fill =colors, cex=0.5)
#add lines
for (i in 1:ncol(Tableplot)) {
  lines(Tableplot[i], type='l', lwd=1.5, col=colors[[i]])
}
dev.off() 
