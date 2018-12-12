args <- commandArgs(TRUE)

library(distory)
library(phangorn)
library(heatmap3)
library(gplots)

trees<-read.tree("gnTrees_collection.tre")
clnum<-length(trees)

geod<-dist.multiPhylo(trees)
wrfd<-wRF.dist(trees)

wrfd<-as.matrix(wrfd)
tagname<-rownames(wrfd)
geod<-as.matrix(geod)
rownames(geod)<-tagname
colnames(geod)<-tagname
geod<-as.dist(geod)

geod.hclust<-hclust(geod, method = "ward.D2")
svg("dendrogram.svg", width=10, height=10)
plot(geod.hclust, cex = 0.6)
dev.off()

gdm<-as.matrix(geod)
svg("heatmap-1.svg", width=12, height=12)
heatmap3(gdm, Rowv = as.dendrogram(geod.hclust), symm = TRUE)
dev.off()

svg("heatmap-2.svg", width=12, height=12)
heatmap.2(gdm, Rowv = as.dendrogram(geod.hclust), symm = TRUE)
dev.off()

dlpre <- as.data.frame(geod.hclust$labels)
dl<- dlpre$`geod.hclust$labels`[order.dendrogram(as.dendrogram(geod.hclust))]
dl <- as.data.frame(dl)
write.table(dl,file = "dend.order", quote = FALSE, col.names = FALSE, row.names = FALSE)

if (clnum > 15){
  clnum<-15
}
for(j in 1:clnum){
cl<-cutree(geod.hclust, k = j)
   if (j < 10){
      j=paste("0",j,sep= "")
      nj=paste("k",j,".cl",sep = "")}
   else
   {nj=paste("k",j,".cl",sep = "")}
write.table(cl, nj, quote = FALSE, col.names = FALSE)
}


