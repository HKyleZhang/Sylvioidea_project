library(ggplot2)
library(tidyr)
logdata<-read.csv("logLkhforR.csv")
a<-logdata[,2:3]
write.csv(a,file = "a.csv")
a<-read.csv("a.csv")
a$X=c("k01","k02","k03","k04","k05","k06","k07","k08","k09","k10","k11","k12","k13","k14")
a<-as.data.frame(a)
names(a)<-c("cluster_num","Actual.data","Average")

logdata<-t(logdata)
b<-logdata[5:104,]
write.csv(b, file = "b.csv")
b<-read.csv("b.csv")
names(b)<-c("X","k01","k02","k03","k04","k05","k06","k07","k08","k09","k10","k11","k12","k13","k14")
num<-seq(1:100)
b$X=num
bb<-b[,2:14]
write.csv(bb, file = "b.csv")

a<-gather(a,group,dk,Actual.data:Average)
b<-gather(b,cluster_num,dk,k01:k14)
write.csv(a,file = "amd.csv")
write.csv(b,file = "bmd.csv")

a<-read.csv("amd.csv")
a<-as.data.frame(a)
b<-read.csv("bmd.csv")
b<-as.data.frame(b)

bp<-ggplot(b,aes(cluster_num,dk)) + geom_boxplot(data=b,aes(cluster_num,dk), color="grey")
res<- bp + geom_point(data=a[1:14,], aes(cluster_num,dk, color="Actual data")) + geom_line(data=a[1:14,], aes(cluster_num,dk, color="Actual data"), group=1) + geom_point(data=a[15:28,], aes(cluster_num,dk, color="Permutation data Mean")) + geom_line(data=a[15:28,], aes(cluster_num,dk, color="Permutation data Mean"), group=1)
ggsave("permutation_analysis.svg", width = 15, height = 15)
