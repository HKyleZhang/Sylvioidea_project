args <- commandArgs(TRUE)

library(treespace)
library(heatmap3)
library(gplots)

trees <- read.tree("gnTrees_collection.tre")
clnum <- length(trees)

L <- as.numeric(args[1])

res <- treespace(trees, method = "treeVec", nf = 2, lambda = L)
kcd <- res$D

kcd.hclust <- hclust(kcd, method = "ward.D2")
svg("dendrogram.svg", width=10, height=10)
plot(kcd.hclust, cex = 0.6)
dev.off()

kcdm <- as.matrix(kcd)
svg("heatmap-1.svg", width=12, height=12)
heatmap3(kcdm, Rowv = as.dendrogram(kcd.hclust), symm = TRUE)
dev.off()

svg("heatmap-2.svg", width=12, height=12)
heatmap.2(kcdm, Rowv = as.dendrogram(kcd.hclust), symm = TRUE)
dev.off()

dlpre <- as.data.frame(kcd.hclust$labels)
dl <- dlpre$`kcd.hclust$labels`[order.dendrogram(as.dendrogram(kcd.hclust))]
dl <- as.data.frame(dl)
write.table(dl,file = "dend.order", quote = FALSE, col.names = FALSE, row.names = FALSE)

if (clnum > 15){
  clnum<-15
}
for(j in 1:clnum){
cl <- cutree(kcd.hclust, k = j)
   if (j < 10){
      j=paste("0",j,sep= "")
      nj=paste("k",j,".cl",sep = "")}
   else
   {nj=paste("k",j,".cl",sep = "")}
write.table(cl, nj, quote = FALSE, col.names = FALSE)
}
