#overview on 6 vs 16 degrees overlap
#use data created in WGCNA script (environment saved)
library(WGCNA)

#data prepearing
networks<- list(expr16C.net, expr6C.net)
moduleLabels<- lapply(networks, function(x) x$colors)
moduleColors<-lapply(networks, function(x) labels2colors(x$colors))
MEs<-lapply(networks, function(x) x$MEs)
MEs<-lapply(MEs, function(x) orderMEs(x, greyName = 'MEO'))
geneTree<-lapply(networks, function(x) x$dendrograms[[1]])

# Isolate the module labels in the order they appear in ordered module eigengenes
ModuleLabels<-lapply(MEs, function(x) substring(names(x), 3))

# Convert the numeric module labels to color labels
Modules<-lapply(ModuleLabels, function(x) labels2colors(as.numeric(x)))

# Numbers of modules for both
nModules<-lapply(Modules, function (x) length(x))

# Initialize tables of p-values and of the corresponding counts
pTable<- matrix(0, nrow = nModules[[1]], ncol = nModules[[2]])
CountTbl<- matrix(0, nrow = nModules[[1]], ncol = nModules[[2]])

# Execute all pairwaise comparisons
for (mod16 in 1:nModules[[1]]) {
  for (mod6 in 1:nModules[[2]]) {
    Members16 = (moduleColors[[1]] == Modules[[1]][mod16])
    Members6 = (moduleColors[[2]] == Modules[[2]][mod6])
    pTable[mod16, mod6] = -log10(fisher.test(Members16, Members6, alternative = "greater")$p.value);
    CountTbl[mod16, mod6] = sum(moduleColors[[1]] == Modules[[1]][mod16] & moduleColors[[2]] == Modules[[2]][mod6])
  }
}


#create color-coded table of the intersection counts
# Truncate p values smaller than 10^{-50} to 10^{-50}
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
               xLabels = paste(" ", Modules[[2]]),
               yLabels = paste(" ", Modules[[1]]),
               colorLabels = TRUE,
               xSymbols = paste("6 ", Modules[[2]], ": ", ModTotal[[2]], sep=""),
               ySymbols = paste("16 ", Modules[[1]], ": ", ModTotal[[1]], sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of 16 set-specific and 6 set-specific modules",
               cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE);
dev.off()

#creating table with only overlaps
ClearTable<- data.frame('NameModule16'= character(), 'NameModule6'=character(), 'N_totgenes16'=integer(), 'N_totgenes6'=integer(), 'N_genestotoverlapped16'= integer(), 'N_genestotoverlapped6'= integer(), 'N_overlapped_genes'=integer(), 'p-value'=integer())
for (mod16 in 1:nModules[[1]]) {
  for (mod6 in 1:nModules[[2]]) {
    Members16 = (moduleColors[[1]] == Modules[[1]][mod16])
    Members6 = (moduleColors[[2]] == Modules[[2]][mod6])
    pTable[mod16, mod6] = -log10(fisher.test(Members16, Members6, alternative = "greater")$p.value);
    if ((sum(moduleColors[[1]] == Modules[[1]][mod16] & moduleColors[[2]] == Modules[[2]][mod6])) >0) {
    lineiwant<-data.frame('NameModule16'= Modules[[1]][mod16], 
                          'NameModule6'= Modules[[2]][mod6], 
                          'N_totgenes16'= length(colnames(expr16C)[moduleColors[[2]]== Modules[[1]][mod16]]), 
                          'N_totgenes6'= length(colnames(expr6C)[moduleColors[[2]]== Modules[[2]][mod6]]) , 
                          'N_genestotoverlapped16'= ModTotal[[1]][mod16], 
                          'N_genestotoverlapped6'= ModTotal[[2]][mod6], 
                          'N_overlapped_genes'=sum(moduleColors[[1]] == Modules[[1]][mod16] & moduleColors[[2]] == Modules[[2]][mod6]), 
                          'p-value'=pTable[mod16,mod6])
    ClearTable<- rbind(ClearTable, lineiwant) 
      }
  }
}


#Bonferroni correction of the pvalue transhold
transhold<- -log10(0.05/6308)

#clean more the table based on transhold
ClearTable2<- ClearTable[ClearTable$p.value > transhold,]

#write.table(ClearTable2, file = '/Volumes/nordborg/pub/forPieter/WGCNA/Results/Table overlaps modules 6 vs 16.txt', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)

#create histograms 
#how many genes are overlapped on the total for that module
hist(ClearTable2$N_overlapped_genes/ClearTable2$N_genes16, breaks = 100, col = 'grey')
hist(ClearTable2$N_overlapped_genes/ClearTable2$N_genes6, breaks = 100, col = 'grey')

ClearTable3<- ClearTable[ClearTable$p.value < transhold,]
length(unique(ClearTable3$NameModule16))
length(unique(ClearTable3$NameModule6))


