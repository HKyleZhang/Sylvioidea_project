out.file <- ""
file.names <- dir(pattern = "logLhcollection")
for(i in 1:length(file.names)){
  file<-read.table(file.names[i])
  out.file <- cbind(out.file, file)
}
write.csv(out.file, file = "logLkhcollection.csv")

da<-read.csv("logLkhcollection.csv")
acdata<-read.table("AutoWZgn_concatenation_logLkhcollection", header = FALSE)
da<-cbind(acdata,da[,3:102])
pernum<-seq(1:100)
names(da)<-c("Actual.data", pernum)
rownum <- nrow(da)

pername<- rep(" ", rownum)
kindex<-seq(1: rownum)
out.file<-cbind(kindex,da$Actual.data,pername,da[,2:101])
names(out.file)<-c("k","Actual data","Permutation data",pernum)
write.csv(out.file, file = "logLkhcollection.csv", row.names = FALSE)

for(k in 1:(rownum - 1)){
  da[k,]<-da[k,]-da[k+1,]
}
da<-rbind(da[1:(rownum - 1),])
dkindex<-seq(1:(rownum - 1))
peravrg<-apply(cbind(da[,2:101]), 1, mean)
persd<-apply(cbind(da[,2:101]), 1, sd)
daforR<-cbind(dkindex,da$Actual.data,peravrg,persd,da[,2:101])
names(daforR)<-c("dk","Actual data","Average","SD",pernum)
write.csv(daforR,file = "logLkhforR.csv", row.names = FALSE)
