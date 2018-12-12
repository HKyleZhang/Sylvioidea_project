#!/bin/bash

#Functions
usage()
{
 echo "Usage: [-a] <Alignment folder> [-L] <A lambda between 0 and 1> [-d] <A different lambda value for cluster assignment> [-m] <cluster method: "ward", "kmeans", "MDSK", "pam", "MDSP", "NMDSK"> [-c] <use codon model> [-k]<method of finding potential number of cluster> [-t] <number of threads> [-s] <seed number for all tree construction>
 [-a] and [-L] are mandatory
 [-c] [-d] [-m] [-t] [-s] are optional"
}

aftertreebuild()
{
 rm -rf ${1}/${2}_rooted
 mkdir ${1}/${2}_rooted
 #mkdir ${1}/treeXtract
 #mkdir ${1}/treeXtract/${2}
 #mkdir ${1}/treeXtract/${2}_rooted

 for file in ${1}/${2}/*treefile
 do
   namepre=$(echo "${file}")
   name=$(basename ${namepre})
   nw_reroot ${file} ZF > ${1}/${2}_rooted/${name}
 done
 
 if [ "${3}" -gt 1 ]; then
    #for file in ${1}/${2}/*treefile
    #do
    #   namepre=$(echo "${file}")
    #   name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "-" -f 2)
    #   seq=$(sed -n "1p" ${file})
    #  
    #   echo "$name $seq" >> ${1}/treeXtract/${2}/gnTrees_collection.tre
    #   echo "$seq" >> ${1}/treeXtract/${2}/gnTrees_noname.tre
    #done

    for file in ${1}/${2}_rooted/*treefile
    do
      namepre=$(echo "${file}")
      if [ "${3}" -gt 5 ]; then
         name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "_" -f 3 | cut -d "-" -f 1)
      else    
         name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "-" -f 2)
      fi 
      seq=$(sed -n "1p" ${file}) 
 
      echo "$name $seq" >> ${1}/gnTrees_collection.tre
   #  echo "$name $seq" >> ${1}/treeXtract/${2}_rooted/gnTrees_collection.tre
   #  echo "$seq" >> ${1}/treeXtract/${2}_rooted/gnTrees_noname.tre
   done
 fi
}

makedatasets()
{
#i=$((RANDOM/37+1))
i=1
#echo -e "The seed for the generation of the 1000-replicate permutation data is: ${3}.\nThe generated datasets(100 datasets in total) start from No.${i} dataset." > impMsg.txt
echo -e "The seed for the generation of the 1000-replicate permutation data is: ${3}." > impMsg.txt

j=$((i+100))
filenum=1
percent=5
rm -rf ${1}/permutation_data
mkdir ${1}/permutation_data  
 
until [ $i -eq $j ]
do
    k=$((i+1))
    startnum=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${i}p")
    end=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${k}p")
    endnum=$((end-1))

    cat ${1}/outfile | sed -n "${startnum},${endnum}p" > ${1}/permutation_data/d${filenum}.phy
    i=$((i+1))
    filenum=$((filenum+1))
done
}

pl_partition()
{
 seqpos1=1
 seqt=0
 for file in ${1}/${2}/*phy
 do 
   name=$(basename "${file}" | cut -d "." -f 1 | cut -d "-" -f 2)
   seq=$(cat ${file} | cut -d " " -f 3)
   seqt=$(echo "$((seqt+seq))")
   echo "${name} = $seqpos1-$seqt;" >> ${1}/Analysis/${3}/partition
   seqpos1=$(echo "$((seqt +1))")
 done
}

partSort()
{
 op=$(echo "${1}_output")
   
 lnum=$(cat partition | wc -l)
 pi=1
 while [ "${pi}" -le "${lnum}" ]
 do
   gnname[$pi]=$(cat partition | sed -n "${pi}p" | cut -d " " -f 1)
   part[$pi]=$(cat partition | sed -n "${pi}p" | cut -d " " -f 3 | cut -d ";" -f 1)
   ((pi++))
 done

 cli=1
 while [ "${cli}" -le "${lnum}" ]
 do
   cl[$cli]=$(cat ${1} | grep "${gnname[$cli]}" | cut -d " " -f 2)
   ((cli++))
 done

 rm -rf ${op}
 i=1
 k=1
 while [ "${i}" -le "${lnum}" ]
 do
   if [ "${i}" -eq "${lnum}" ]; then
      j=${i}
   else
      j=$((i+1))
   fi

   if [ "${cl[$i]}" != "NA" ]; then
      partcon=$(echo "${part[$i]}")
      while [ "${j}" -le "${lnum}" ]
      do
        if [ "${cl[$j]}" != "NA" ]; then
           if [ "${cl[$j]}" -eq "${cl[$i]}" ]; then
              if [ "${i}" -lt "${lnum}" ]; then
                 partcon=$(echo "${partcon},${part[$j]}")
                 cl[$j]="NA"
              fi
           fi
        fi

        ((j++))
      done  
      if [ "${k}" -lt 10  ]; then
         echo "cluster0${k} = ${partcon}" >> ${op}
      else
         echo "cluster${k} = ${partcon}" >> ${op}
      fi
      ((k++))
   fi
   ((i++))
 done
}

logLXtract()
{
  clfiup=$(ls ${2}/*iqtree | wc -l)
  k=1
  while [ "${k}" -le "${clfiup}" ]
  do
    if [ "${k}" -lt 10 ]; then kmd=$(echo "0${k}"); else kmd=${k}; fi
    logLh=$(cat ${2}/*cluster${kmd}*iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
    logLhmd=$(echo "${logLh} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
    ttllogLhmd=$((ttllogLhmd + logLhmd))
    ((k++))    
  done
  seqname=$(echo "${2}" | rev | cut -d "/" -f 2 | rev)
  echo -e "${seqname} k=${clfiup} -loglikelihood*10^4:\t${ttllogLhmd}" >> ${1}/${seqname}_logLhcollection
}

MinFinder()
{
 fnum=$(ls ${1} | wc -l)
 alnname=$(echo "${1}" | rev | cut -d "/" -f 1 | rev)
 bsreppre=$(cat ${1}/${alnname}-cluster01.kcd | wc -l)
 bsrep=$((bsreppre - 1))
 i=1
 while [ "${i}" -le "${bsrep}" ]
 do
   j=1
   while [ "${j}" -le "${fnum}" ]
   do
      if [ "${j}" -lt 10 ]; then jmd=$(echo "cluster0${j}"); else jmd=$(echo "cluster${j}"); fi       
      kcdvpre=$(cat ${1}/${alnname}-${jmd}.kcd | sed -n "${i}p")
      kcdv[$j]=$(echo "${kcdvpre} * 100000000" | bc | cut -d "." -f 1)
      ((j++))
   done

   j=1
   kcdvmin=${kcdv[1]}
   min=1
   while [ "${j}" -le "${fnum}" ]
   do
     if [ "${kcdv[$j]}" -lt "${kcdvmin}" ]; then
        kcdvmin=${kcdv[$j]}
        min=${j}
     fi
     ((j++))
   done
       
   echo "${min}" >> ${2}/${alnname}-cluster.assign
 
   ((i++))
 done
}


#R implementation
get_cluster='if (!require(treespace)) install.packages("treespace")
if (!require(heatmap3)) install.packages("heatmap3")
if (!require(gplots)) install.packages("gplots")
if (!require(cluster)) install.packages("cluster")
if (!require(smacof)) install.packages("smacof")

trees <- read.tree("gnTrees_collection.tre")
clnum <- length(trees)

L <- args1

res <- treespace(trees, method = "treeVec", nf = 2, lambda = L, return.tree.vectors = TRUE)

if ("args2" == "ward") {
   kcd <- res$D
   
   #Multi-Dimensional Scaling for convenience of visualzing the results
   mdres2D <- smacofSym(kcd, ndim = 2, type = "ratio")
   treMDS2D <- as.data.frame(mdres2D$conf)
   write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
   write.table(mdres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    
   mdres3D <- smacofSym(kcd, ndim = 3, type = "ratio")
   treMDS3D <- as.data.frame(mdres3D$conf)
   write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
   write.table(mdres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)

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

   if (clnum > 15){clnum<-15}

   for(j in 1:clnum){
      cl <- cutree(kcd.hclust, k = j)
       if (j < 10){
          j=paste("0",j,sep= "")
          nj=paste("k",j,".cl",sep = "")
        } else {
          nj=paste("k",j,".cl",sep = "")
        }
      write.table(cl, nj, quote = FALSE, col.names = FALSE)
    }
}


# K-means on tree vectors
if ("args2" == "kmeans") {
   treMDS <- res$pco$li
   write.table(treMDS, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
   
   treVec <- res$vectors

   if (clnum > 15){
   	  clnum <- 15
    } else {
      clnum <- clnum - 1
    }

   for (j in 1: clnum){
       kc.km <- kmeans(treVec, centers = j, nstart = 50)
       if (j < 10){
          j=paste("0",j,sep= "")
          nj=paste("k",j,".cl",sep = "")
        } else {
          nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
    }
}

# Multi-Dimensional scaling down to 2D and 3D, K-means on 3D MDS
if ("args2" == "MDSK") {
	kcd <- res$D
    mdres2D <- smacofSym(kcd, ndim = 2, type = "ratio")
    treMDS2D <- as.data.frame(mdres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    
    mdres3D <- smacofSym(kcd, ndim = 3, type = "ratio")
    treMDS3D <- as.data.frame(mdres3D$conf)
    write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)

    if (clnum > 15){
       clnum <- 15
    } else {
       clnum <- clnum - 1
    }

    for (j in 1: clnum){
        kc.km <- kmeans(treMDS3D, centers = j, nstart = 50)
        if (j < 10){
           j=paste("0",j,sep= "")
           nj=paste("k",j,".cl",sep = "")
        } else {
           nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
    }
}

# PAM on the distance/dissimilarities matrix
if ("args2" == "pam") {
   kcd <- res$D

   if (clnum > 15){
   	  clnum <- 15
    } else {
      clnum <- clnum - 1
    }

   for (j in 1: clnum){
       kc.pam <- pam(kcd, diss = TRUE, k = j, cluster.only = TRUE)
       if (j < 10){
          j=paste("0",j,sep= "")
          nj=paste("k",j,".cl",sep = "")
        } else {
          nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.pam, nj, quote = FALSE, col.names = FALSE)
    }
}

# Multi-Dimensional scaling down to 2D and 3D, PAM on 3D MDS
if ("args2" == "MDSP") {
	kcd <- res$D
    mdres2D <- smacofSym(kcd, ndim = 2, type = "ratio")
    treMDS2D <- as.data.frame(mdres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    
    mdres3D <- smacofSym(kcd, ndim = 3, type = "ratio")
    treMDS3D <- as.data.frame(mdres3D$conf)
    write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)

    if (clnum > 15){
       clnum <- 15
    } else {
       clnum <- clnum - 1
    }

    for (j in 1: clnum){
        kc.pam <- pam(treMDS3D, diss = FALSE, k = j, cluster.only = TRUE)
        if (j < 10){
           j=paste("0",j,sep= "")
           nj=paste("k",j,".cl",sep = "")
        } else {
           nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.pam, nj, quote = FALSE, col.names = FALSE)
    }
}

# Non-metric Multi-Dimensional scaling down to 2D and K-means
if ("args2" == "NMDSK") {
	kcd <- res$D
    mdres2D <- smacofSym(kcd, ndim = 2, type = "ordinal")
    treMDS2D <- as.data.frame(mdres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)

    if (clnum > 15){
       clnum <- 15
    } else {
       clnum <- clnum - 1
    }

    for (j in 1: clnum){
        kc.km <- kmeans(treMDS2D, centers = j, nstart = 50)
        if (j < 10){
           j=paste("0",j,sep= "")
           nj=paste("k",j,".cl",sep = "")
        } else {
           nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
    }
}'

logLkhforR_make='out.file <- ""
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
write.csv(daforR,file = "logLkhforR.csv", row.names = FALSE)'

graphs='if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")

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
ggsave("permutation_analysis.svg", width = 15, height = 15)'

find_potent_knum='logd <- read.csv("logLkhforR.csv", header = TRUE)
logd <- abs(logd[,2]-logd[,3])
logdmin <- min(logd, na.rm = TRUE)
min.ind <- which(logd == logdmin)

logper <- seq(1:14)
for (i in 1:13) {
  logper[i+1] <- logd[i+1]/logd[1]
}


res <- array("NA", 14)
j <- 1
threshd <- 0.5
for (i in 1:12){
  if (i >= 2){
    if (logper[i] <= threshd){
      a <- logd[i-1] - logd[i]
      b <- logd[i] - logd[i+1]
      if (b <= a) {
        if (logper[i+1] <= threshd){
        res[j] <- i
        j <- j+1
        }
      }
    }
  }
  else {
    if (logd[i] == logdmin){
      res[j] <- i
      j <- j+1
    }
  }
}
res <- res[res != "NA"]
res <- as.data.frame(res)
write.table(res, file = "potent_knum", quote = FALSE, col.names = FALSE, row.names = FALSE)'

KCmetric='if (!require(treespace)) install.packages("treespace")
if (!require(foreach)) install.packages("foreach")
if (!require(doParallel)) install.packages("doParallel")

no_cores <- detectCores()
registerDoParallel(no_cores)

L <- args1

files <- list.files(pattern = "*tre")
foreach(i=files) %dopar% {
  print(i)
  name <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  name <- paste(name, ".kcd", sep = "")
  trees <- read.tree(i)
  res <- multiDist(trees, lambda = L)
  res <- as.matrix(res)
  kcd <- round(res[2:length(res[,1]),1],8)
  write.table(kcd, file = name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

stopImplicitCluster()'

clusterAssign='pknum <- as.numeric(args1)

if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(svglite)) install.packages("svglite")

data <- read.table("Cluster_Assignment.summary", header = FALSE)
datash <- read.table("Cluster_Assignment_shortdata.summary", header = FALSE)

names(data) <- c("order", "gene","cluster","percentage")
data$percentage <- data$percentage / 10
data$order <- factor(data$order)

if (pknum == 4){
  names(datash) <- c("order", "gene","Cluster1","Cluster2","Cluster3","Cluster4")
  res<-ggplot() + geom_bar(data=data, aes(y = percentage,x = order,fill = cluster), stat="identity", position = "stack", width = 0.9) + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, angle = 90), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm")) + scale_x_discrete(position = "bottom", labels = datash$gene) + scale_fill_manual(values = c("coral2", "aquamarine4", "yellow3", "deepskyblue2"))
}
if (pknum == 5){
  names(datash) <- c("order", "gene","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")
  res<-ggplot() + geom_bar(data=data, aes(y = percentage,x = order,fill = cluster), stat="identity", position = "stack", width = 0.9) + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, angle = 90), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm")) + scale_x_discrete(position = "bottom", labels = datash$gene) + scale_fill_manual(values = c("coral2", "deepskyblue2", "yellow3", "aquamarine4", "blueviolet"))
}
if (pknum == 6){
  names(datash) <- c("order", "gene","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6")
  res<-ggplot() + geom_bar(data=data, aes(y = percentage,x = order,fill = cluster), stat="identity", position = "stack", width = 0.9) + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, angle = 90), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm")) + scale_x_discrete(position = "bottom", labels = datash$gene) + scale_fill_manual(values = c("coral2", "deepskyblue2", "yellow3", "aquamarine4", "blueviolet", "magenta2"))
}
if (pknum == 7){
  names(datash) <- c("order", "gene","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7")
  res<-ggplot() + geom_bar(data=data, aes(y = percentage,x = order,fill = cluster), stat="identity", position = "stack", width = 0.9) + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, angle = 90), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm")) + scale_x_discrete(position = "bottom", labels = datash$gene) + scale_fill_manual(values = c("coral2", "aquamarine4", "blueviolet", "magenta2", "deepskyblue2", "yellow3", "gray60"))
}

ggsave("Cluster_Assignment.svg", width = 15, height = 5)'

editor='if (!require(foreach)) install.packages("foreach")
if (!require(doParallel)) install.packages("doParallel")

no_cores <- detectCores()
registerDoParallel(no_cores)

files <- list.files(pattern = "*model")
foreach(i=files) %dopar% {
    name <- strsplit(i, ".", fixed = TRUE)[[1]][1]
    name <- paste(name, ".txt", sep = "")
  dd <- read.table(i)
  optmodel <- dd[1,2]
  for (k in 1:nrow(dd)){
    dd[k,2] <- round(((dd[k,2] - optmodel) ^ 2), 3)
  }
  write.table(dd[order(dd$V1),], file = name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}'

find_model='if (!require(dplyr)) install.packages("dplyr")
dd <- read.table("model.summary")
dd$sum_sq <- apply(dd[,2:ncol(dd)], 1, sum)
dd <- arrange(dd, sum_sq)
write.table(dd, "bestmodel.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
'



#--------------------------------------------------------------------------------#
#Loading the flags
while getopts "ha:c:u:L:d:m:k:t:s:" opt
do
  case $opt in
               h) usage && exit;;
               a) alnf=$(echo "${OPTARG}" | cut -d "/" -f 1);;
               c) cod="CODON";;
               L) lambda=${OPTARG};;
               d) lambda2=${OPTARG};;
               m) clustermethod=${OPTARG};;
               k) optik=${OPTARG};;
               t) thr=${OPTARG};;
               s) iqseed=${OPTARG};;
   esac
done

if [ $OPTIND -eq 1 ]; then
   usage
   exit
fi


curdir=$(pwd)
resize -s 39 80 > /dev/null
if [ -z "${thr}" ]; then thr=$(nproc --all); fi
if [ -z "${lambda}" ]; then 
   if [ -z "${lambda2}" ]; then
      echo -e "\nERROR: Missing lambda value!" && exit
   fi
else
   if [ -z "${lambda2}" ]; then
      lambda2=${lambda}
   fi
fi

if [ -z "${clustermethod}" ]; then clustermethod="ward"; fi	

if [ -z "${iqseed}" ]; then iqseed=${RANDOM}; fi


if [ -z "${optik}" ]; then
   case ${clustermethod} in
   	                       ward) optik="gori";;
                         kmeans) optik="clustree";;
                           MDSK) optik="clustree";;
                            pam) optik="gori";;
                           MDSP) optik="clustree";;
                          NMDSK) optik="clustree";;
   esac
fi

if [ -n "${cod}" ]; then
   mk=92
   mf="-st CODON"
   codonmodel="Yes"
else
  mk=134
  mf=""
  codonmodel="No"
fi

echo '
          _____                                    
         / ____\          /|           _         _      /|
        / /    \|        | |          /_|       /_|    | |
        \ \__  __      __| |__      __ _   ___   _   __| |  ___    ____
         \___ \\ \    / /| |\ \    / /| | / _ \ | | / _  | / _ \  / _  \  
      _      \ \\ \  / / | | \ \  / / | || / \ || || / | || |_| || / | | 
      \\_____/ / \ \/ /  | |  \ \/ /  | || \_/ || || \_| |\  __ || \_| |  
       \______/   \  /   |_|   \__/   |_| \___/ |_| \_____\\_\    \_____\
==================/ /===========================================================
=================/ /==================================== U-l-t-i-m-a-t-e =======
                /_/                                                         '
read -p "
 Files checklist: 1.Alignment folder [-a]

    Dependencies: 1.IQ-Tree               4.fseqboot
                  2.Newick Utilities      5.Rstudio
                  3.AMAS                  

      R packages: 1.treespace             7.dplyr
                  2.heatmap3              8.foreach
                  3.gplots                9.doParallel
                  4.ggplot2               10.clustree
                  5.tidyr                 11.svglite
                  6.smacof
   *************************************************************
     SETTINGS:

                alignment folder: ${alnf}
                     codon model: ${codonmodel}
                          lambda: ${lambda}
                         lambda2: ${lambda2}
                  cluster method: ${clustermethod}
                  find optimal k: ${optik}
                         threads: ${thr}
                            seed: ${iqseed}
 
   *************************************************************

Continue?(Y/n): " answer

if [ "${answer}" == "Y" -o "${answer}" == "y" ]; then
   if [ -z "${alnf}" ] || [ ! -d "${curdir}/${alnf}" ]; then echo -e "\nMsg: Some files are missing.\n" && usage && exit; fi  #Check if all the files exist 

echo '
          _____                                    
         / ____\          /|           _         _      /|
        / /    \|        | |          /_|       /_|    | |
        \ \__  __      __| |__      __ _   ___   _   __| |  ___    ____
         \___ \\ \    / /| |\ \    / /| | / _ \ | | / _  | / _ \  / _  \  
      _      \ \\ \  / / | | \ \  / / | || / \ || || / | || |_| || / | | 
      \\_____/ / \ \/ /  | |  \ \/ /  | || \_/ || || \_| |\  __ || \_| |  
       \______/   \  /   |_|   \__/   |_| \___/ |_| \_____\\_\    \_____\
==================/ /===========================================================
=================/ /==================================== U-l-t-i-m-a-t-e =======
                /_/                                                         ' > ${curdir}/parameter_input

   if [ "${clustermethod}" == "MDSK" ]; then
      echo "  SETTINGS: 
      
                alignment folder: ${alnf}
                     codon model: ${codonmodel}
                          lambda: ${lambda}
                         lambda2: ${lambda2}
                  cluster method: MDS(Multi-Dimensional Scaling) -> K-means
                  find optimal k: ${optik}
                            seed: ${iqseed}

  Time stamp: $(date)" >> ${curdir}/parameter_input
   elif [ "${clustermethod}" == "MDSP" ]; then
   	  echo "  SETTINGS: 

                alignment folder: ${alnf}
                     codon model: ${codonmodel}
                          lambda: ${lambda}
                         lambda2: ${lambda2}
                  cluster method: MDS(Multi-Dimensional Scaling) -> PAM
                  find optimal k: ${optik}
                            seed: ${iqseed}

  Time stamp: $(date)" >> ${curdir}/parameter_input
   elif [ "${clustermethod}" == "NMDSK" ]; then
   	  echo "  SETTINGS: 
                
                alignment folder: ${alnf}
                     codon model: ${codonmodel}
                          lambda: ${lambda}
                         lambda2: ${lambda2}
                  cluster method: NMDS(non-metric Multi-Dimensional Scaling) -> K-means
                  find optimal k: ${optik}
                            seed: ${iqseed}

  Time stamp: $(date)" >> ${curdir}/parameter_input
   else
	 echo "  SETTINGS: 
                
                alignment folder: ${alnf}
                     codon model: ${codonmodel}
                          lambda: ${lambda}
                         lambda2: ${lambda2}
                  cluster method: ${clustermethod}
                  find optimal k: ${optik}
                            seed: ${iqseed}

  Time stamp: $(date)" >> ${curdir}/parameter_input
   fi

#read -p $'\nSpecify the directory of AMAS.py (Default is $HOME/Software/amas) ' amasdir
   amasdir="$HOME/Software/amas"

tstamp=$(date)
echo "
Analysis started at ${tstamp}
--------------------------------------------------------------------------------"

partindex=1
#1.Building gene trees
#*************************************************************
echo "
Part ${partindex}. Gene trees construction"
SECONDS=0

fdname="iqMF-autoWZgn_1000UFB"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia ${curdir}/${fdname}/bs-file

bsset="-bb 1000 -bnni -wbtl" #Perform Ultra Fast Bootstrapping with 1000 replicates.
cd ${curdir}/${alnf}/
parallel --no-notice -j ${thr} "iqtree -s {} ${mf} -AICc -ninit 200 -ntop 50 -nt 1 -seed ${iqseed} ${bsset} -keep-ident -quiet" ::: *.phy

#File sorting
mv ${curdir}/${alnf}/*ufboot ${curdir}/${fdname}/bs-file
mv ${curdir}/${alnf}/*iqtree ${curdir}/${fdname}/report/
mv ${curdir}/${alnf}/*treefile ${curdir}/${fdname}/
mv ${curdir}/${alnf}/*log ${curdir}/${alnf}/*gz ${curdir}/${alnf}/*bionj ${curdir}/${alnf}/*mldist ${curdir}/${alnf}/*contree ${curdir}/${alnf}/*splits.nex ${curdir}/${fdname}/trivia
if [ -n "${cod}" ];then mv ${curdir}/${alnf}/*parstree ${curdir}/${fdname}/trivia; fi

#Extract the model information
for i in ${curdir}/${fdname}/report/*iqtree
do
  name=$(basename ${i} | cut -d "." -f 1-2)
  model=$(cat ${i} | grep "Best-fit model according to AICc:" | cut -d " " -f 6)
  echo -e "${name}\t${model}" >> ${curdir}/${alnf}/iq-modelset
done

#Find a single model for permutation analysis
rm -rf ${curdir}/report-temp
mkdir ${curdir}/report-temp
for i in ${curdir}/${fdname}/report/*iqtree
do
  name=$(basename ${i} | cut -d "." -f 1)
  cat ${i} | grep -A${mk} "List of models sorted by AICc scores:" > ${curdir}/report-temp/${name}
done

rm -rf ${curdir}/report-temp2
mkdir ${curdir}/report-temp2
for i in ${curdir}/report-temp/*
do
  name=$(basename ${i})
  lnum=$(cat ${i} | wc -l)
  sed -n "4,${lnum}p" ${i} > ${curdir}/report-temp2/${name}
done
rm -rf ${curdir}/report-temp

rm -rf ${curdir}/report-model
mkdir ${curdir}/report-model
for i in ${curdir}/report-temp2/*
do
  name=$(basename ${i})
  cat ${i} | cut -c 1-17,30-39 > ${curdir}/report-model/${name}.model
done
rm -rf ${curdir}/report-temp2
cd ${curdir}/report-model/
Rscript <(echo "${editor}") 2>&1 >/dev/null
rm -rf ${curdir}/report-model/*model

rm -rf ${curdir}/report-temp
mkdir ${curdir}/report-temp
for i in ${curdir}/report-model/*
do 
   if [ ! -e ${curdir}/report-temp/00A.txt ]; then
      cat ${i} | cut -d " " -f 1 > ${curdir}/report-temp/00A.txt
   fi
   name=$(basename ${i})
   cat ${i} | cut -d " " -f 2 > ${curdir}/report-temp/${name}
done
rm -rf ${curdir}/report-model

cd ${curdir}/report-temp/
paste -d " " * > ${curdir}/model.summary
rm -rf ${curdir}/report-temp

cd ${curdir}
Rscript <(echo "${find_model}") 2>&1 >/dev/null
rm -rf ${curdir}/model.summary

kmodel=$(cat ${curdir}/bestmodel.txt | sed -n "1p" | cut -d " " -f 1)

#Root the tree and extract the rooted trees into a file.
stepind=1
aftertreebuild ${curdir} ${fdname} ${stepind}

cp -r ${curdir}/${fdname}/bs-file ${curdir}/

rm -rf ${curdir}/Trees
mkdir ${curdir}/Trees
mv ${curdir}/${fdname}/ ${curdir}/${fdname}_rooted/ ${curdir}/Trees/

#Reroot the bootstrapped trees
rm -rf ${curdir}/bs-file_rooted
mkdir ${curdir}/bs-file_rooted
for ufbfile in ${curdir}/bs-file/*
do
  namepre=$(echo "${ufbfile}")
  name=$(basename ${namepre})
  nw_reroot ${ufbfile} ZF > ${curdir}/bs-file_rooted/${name}
done

rm -rf ${curdir}/bs-file/
cd ${curdir}

((partindex++))
duration=$SECONDS
tstamp=$(date)
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${tstamp}

--------------------------------------------------------------------------------"
#*************************************************************

#2.Generate permutation data/concatenation sequence
#*************************************************************
if [ "${optik}" == "gori" ]; then
   echo "
Part ${partindex}. Generate permutation data"
   SECONDS=0

   python ${amasdir}/AMAS.py concat -i ${alnf}/*phy -f phylip -d dna -u phylip -t ${curdir}/AutoWZgn_concatenation.phy > /dev/null
   rm -rf ${curdir}/partitions.txt
   sed -i "s/ZF/ZF    /" AutoWZgn_concatenation.phy
   sed -i "s/MW/MW    /" AutoWZgn_concatenation.phy
   sed -i "s/CW/CW    /" AutoWZgn_concatenation.phy
   sed -i "s/BC/BC    /" AutoWZgn_concatenation.phy
   sed -i "s/GRW/GRW    /" AutoWZgn_concatenation.phy
   sed -i "s/CRW/CRW    /" AutoWZgn_concatenation.phy

   seednum=2
   until [ $((seednum % 2)) -gt 0 ]
   do
     seednum=${RANDOM}
   done
   fseqboot -sequence ${curdir}/AutoWZgn_concatenation.phy -outfile outfile -test o -seqtype d -rewriteformat p -reps 101 -seed ${seednum} > /dev/null

   taxanum=6
   makedatasets ${curdir} ${taxanum} ${seednum}
   rm -rf ${curdir}/outfile
else
   echo "
Part ${partindex}. Generate concatenation sequence"
   SECONDS=0

   python ${amasdir}/AMAS.py concat -i ${alnf}/*phy -f phylip -d dna -u phylip -t ${curdir}/AutoWZgn_concatenation.phy > /dev/null
   rm -rf ${curdir}/partitions.txt
   sed -i "s/ZF/ZF    /" AutoWZgn_concatenation.phy
   sed -i "s/MW/MW    /" AutoWZgn_concatenation.phy
   sed -i "s/CW/CW    /" AutoWZgn_concatenation.phy
   sed -i "s/BC/BC    /" AutoWZgn_concatenation.phy
   sed -i "s/GRW/GRW    /" AutoWZgn_concatenation.phy
   sed -i "s/CRW/CRW    /" AutoWZgn_concatenation.phy
fi

((partindex++))
duration=$SECONDS
tstamp=$(date)
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${tstamp}

--------------------------------------------------------------------------------"
#*************************************************************

#3.Prepare files
#*************************************************************
if [ "${optik}" == "gori" ]; then
    echo "
Part ${partindex}. File preparation for permutation analysis"
else
	echo "
Part ${partindex}. File preparation for clustering analysis"
fi

SECONDS=0

#Generate gene trees collection
fdname="iqMF-autoWZgn_br"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia

bsset="" #Don't need nodes support method.
cd ${curdir}/${alnf}/
parallel --no-notice -j ${thr} --colsep '\t' "iqtree -s {1} ${mf} -m {2} -ninit 200 -ntop 50 -nt 1 -seed ${iqseed} ${bsset} -keep-ident -quiet" :::: iq-modelset
mv iq-modelset ${curdir}/

#File sorting
mv ${curdir}/${alnf}/*iqtree ${curdir}/${fdname}/report/
mv ${curdir}/${alnf}/*treefile ${curdir}/${fdname}/
mv ${curdir}/${alnf}/*log ${curdir}/${alnf}/*gz ${curdir}/${alnf}/*bionj ${curdir}/${alnf}/*mldist ${curdir}/${fdname}/trivia
if [ -n "${cod}" ];then mv ${curdir}/${alnf}/*parstree ${curdir}/${fdname}/trivia; fi

##Root the tree and extract the rooted trees into a file.
stepind=5
aftertreebuild ${curdir} ${fdname} ${stepind}

for trefile in ${curdir}/${fdname}_rooted/*treefile
do
  name=$(basename ${trefile} | cut -d "." -f 1-2)  
  cat ${curdir}/bs-file_rooted/${name}.ufboot ${trefile} > ${curdir}/bs-file_rooted/${name}.ufboot-temp
  rm ${curdir}/bs-file_rooted/${name}.ufboot
  mv ${curdir}/bs-file_rooted/${name}.ufboot-temp ${curdir}/bs-file_rooted/${name}.ufboot
done

mv ${curdir}/${fdname}/ ${curdir}/${fdname}_rooted/ ${curdir}/Trees/

#Get clusters information
rm -rf ${curdir}/Analysis
mkdir ${curdir}/Analysis
mv ${curdir}/gnTrees_collection.tre ${curdir}/Analysis/
if [ -e ${curdir}/impMsg.txt ]; then mv ${curdir}/impMsg.txt ${curdir}/Analysis/; fi

cd ${curdir}/Analysis
get_cluster_temp=$(echo "${get_cluster}" | sed "s/args1/${lambda}/" | sed "s/args2/${clustermethod}/")
Rscript <(echo "${get_cluster_temp}") 2>&1 >/dev/null

#Summarize the cluster information for each gene in each number of cluster
kfnum=$(ls ${curdir}/Analysis/*cl | wc -l)
i=1
header=""
while [ "${i}" -le "${kfnum}" ]
do
  header=$(echo "${header} k${i}")
  ((i++))
done
echo "genes${header}" > ${curdir}/Analysis/cluster.summary

while read line
do
  kf=2
  ksort=""
  while [ "${kf}" -le "${kfnum}" ]
  do
  	if [ "${kf}" -lt 10 ]; then kfmd=$(echo "k0${kf}.cl"); else kfmd=$(echo "k${kf}.cl"); fi
    name=$(echo "${line}" | cut -d " " -f 1)
    ki=$(cat ${curdir}/Analysis/${kfmd} | grep "${name}" | cut -d " " -f 2)
    ksort=$(echo "${ksort} ${ki}")
    ((kf++))
  done
  echo "${line}${ksort}" >> ${curdir}/Analysis/cluster.summary
done < ${curdir}/Analysis/k01.cl

#Generate order file to guide the summary of cluster assignment
if [ -e ${curdir}/Analysis/dend.order ]; then
	sexyn="n"
    lnum=$(cat ${curdir}/Analysis/dend.order | wc -l)
    i=1
    while [ "${i}" -le "${lnum}" ]
    do
      name=$(cat ${curdir}/Analysis/dend.order | sed -n "${i}p")
      nsmark=$(echo "${name}" | cut -c 1-2)
      if [ "${nsmark}" == "ns" ]; then 
         prefix=$(echo "${name}" | cut -d "_" -f 1)
         sexyn="y"
         echo "${name}" >> ${curdir}/Analysis/Neo.sexpre
      else
         echo "${name}" >> ${curdir}/Analysis/A.auto
      fi
      ((i++))
    done
    rm -rf dend.order

    if [ "${sexyn}" == "y" ]; then
       cat ${curdir}/Analysis/Neo.sexpre | sort > ${curdir}/Analysis/Neo.sex
       rm -rf ${curdir}/Analysis/Neo.sexpre
       cat ${curdir}/Analysis/A.auto ${curdir}/Analysis/Neo.sex > ${curdir}/Analysis/order
    else
       mv ${curdir}/Analysis/A.auto ${curdir}/Analysis/order
    fi
    rm -rf ${curdir}/Analysis/A.auto ${curdir}/Analysis/Neo.sex

    if [ ! -d clinfo ]; then mkdir clinfo; fi
    mv *cl order clinfo/
fi

if [ ! -d clinfo ]; then mkdir clinfo && mv *cl clinfo/; fi

if [ "${optik}" == "gori" ]; then
   rm -rf permutationGTRI_analysis
   mkdir permutationGTRI_analysis
   mv clinfo/ permutationGTRI_analysis/

   #Get partition file
   pl_partition ${curdir} ${alnf} permutationGTRI_analysis

   #Manoevour the sequence files
   cp -p ${curdir}/AutoWZgn_concatenation.phy ${curdir}/Analysis/permutationGTRI_analysis/
   mv ${curdir}/permutation_data/*phy ${curdir}/Analysis/permutationGTRI_analysis/
   rmdir ${curdir}/permutation_data
else
  rm -rf optimal_cluster
  mkdir optimal_cluster

  #Get partition file
  pl_partition ${curdir} ${alnf} optimal_cluster
  
  #Manoevour the files
	mv ${curdir}/AutoWZgn_concatenation.phy ${curdir}/Analysis/clinfo/ ${curdir}/Analysis/cluster.summary ${curdir}/Analysis/optimal_cluster/
fi

((partindex++))
duration=$SECONDS
tstamp=$(date)
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${tstamp}

--------------------------------------------------------------------------------"
#*************************************************************

#4.Permutation analysis/Clustering analysis
#*************************************************************
if [ "${optik}" == "gori" ]; then
   echo "
Part ${partindex}. Permutation analysis"
   SECONDS=0

   cd ${curdir}/Analysis/permutationGTRI_analysis/
   predir=$(pwd)

   #Create partition-sorted file
   clfi=0
   rm -rf cl_sort
   mkdir cl_sort
   for file in clinfo/*cl
   do  
     partSort ${file}
     filename=$(basename ${file})
     ((clfi++))
     clf[$clfi]=$(echo "${filename}_output")
     mv ./clinfo/*cl_output ./cl_sort/
   done

   #Execute the sequence splitting
   ni=0
   for file in *phy
   do       
     ((ni++))
     name[$ni]=$(echo "${file}" | cut -d "." -f 1)
     mkdir ${name[$ni]}

     exei=1   
     while [ "${exei}" -le "${clfi}" ]     
     do
       python ${amasdir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ./cl_sort/${clf[$exei]} -u phylip > /dev/null
       kclnum=$(echo "${clf[$exei]}" | cut -d "." -f 1)
       mkdir ./${name[$ni]}/${kclnum}
       mv ./*cluster* ./${name[$ni]}/${kclnum}/
       ((exei++))
     done
     rm -rf ${file}
   done

   #read -p $'\nContinue to Tree Construction? (y/n): ' totree
   totree="y"
   if [ "${totree}" = "y" ] || [ "${totree}" = "Y" ]; then

      tri=1
      while [ "${tri}" -le "${ni}" ]
      do
        trj=1
        while [ "${trj}" -le "${clfi}" ] 
        do    
          kclnum=$(echo "${clf[$trj]}" | cut -d "." -f 1)
          cd ${predir}/${name[$tri]}/${kclnum}/
          clnum=1
          for seqfile in *
          do         
            directory=$(pwd)
            echo "${directory}/${seqfile}" >> ${curdir}/permtdata
            echo "${directory}" >> ${curdir}/permtdir-temp
            done
            ((trj++))
          done
         ((tri++))
       done
   
      cd ${curdir}
      parallel --no-notice -j ${thr} "iqtree -s {} ${mf} -m ${kmodel} -ninit 200 -ntop 50 -nt 1 -seed ${iqseed} -keep-ident -quiet" :::: permtdata
      rm -rf permtdata
   
      cat ${curdir}/permtdir-temp | sort -n | uniq > ${curdir}/permtdir
      rm -rf ${curdir}/permtdir-temp

      export -f logLXtract
      parallel --no-notice -j ${thr} logLXtract ::: ${predir} :::: permtdir > /dev/null
      rm -rf permtdir
      if [ ! -d "${predir}/logLikelihood" ]; then mkdir ${predir}/logLikelihood; fi
      mv ${predir}/*logLhcollection ${predir}/logLikelihood/ 
      rm -rf ${predir}/d*

      #Extract log likelihood values from numerous files:
      mkdir ${predir}/value
      cd ${predir}/logLikelihood
      for file in *
      do
        name=$(echo "${file}")
        cat ${file} | cut -d $'\t' -f 2 > ${predir}/value/${name}
      done

      mv ${predir}/value/AutoWZgn_concatenation_logLhcollection ${predir}/value/AutoWZgn_concatenation_logLkhcollection

      #Generate logLkhforR.csv file
      cd ${predir}/value
      Rscript <(echo "${logLkhforR_make}") 2>&1 >/dev/null

      #File sorting
      mv ${predir}/value ${predir}/logLikelihood
      mv ${predir}/logLikelihood/value/logLkhcollection.csv ${predir}/logLikelihood/value/logLkhforR.csv ${predir}/
      rm -rf ${predir}/logLikelihood

      #Plot graph
      cd ${predir}
      Rscript <(echo "${graphs}") 2>&1 >/dev/null
      cd ${curdir}/Analysis/
      Rscript -e 'library("clustree"); dd <- read.table("cluster.summary", header = TRUE); p <- clustree(dd, prefix = "k"); ggsave("cluster_summary.svg", width = 12, height = 12)' 2>&1 >/dev/null
 
      #File sorting
      rm -rf ${predir}/a.csv ${predir}/amd.csv ${predir}/b.csv ${predir}/bmd.csv ${predir}/cl_sort
   fi
else
     echo "
Part ${partindex}. Find optimal number of cluster"
   SECONDS=0
   cd ${curdir}/Analysis/optimal_cluster/
   predir=$(pwd)
   Rscript -e 'library("clustree"); dd <- read.table("cluster.summary", header = TRUE); p <- clustree(dd, prefix = "k"); ggsave("cluster_summary.svg", width = 12, height = 12)' 2>&1 >/dev/null
   
   if [ "${optik}" == "clustree" ]; then	
      #Visualize the cluster
      sc=1
      while [ "${sc}" -le 1 ]
      do
     	 if [ -e ${curdir}/Analysis/cluster.visualization ]; then rm -rf ${curdir}/Analysis/cluster.visualization; fi

         read -p $'\nCheck the cluster_summary.svg and specify the number of cluser to visualize: ' vk
         if [ "${vk}" -lt 10 ]; then vkmd=$(echo "k0${vk}.cl"); else vkmd=$(echo "k${vk}/cl"); fi
         echo "genes M1 M2 cluster" > ${curdir}/Analysis/cluster.visualization
         while read line
         do
        	name=$(echo "${line}" | cut -d " " -f 1)
        	k=$(cat ${predir}/clinfo/${vkmd} | grep "${name}" | cut -d " " -f 2)
         	echo "${line} ${k}" >> ${curdir}/Analysis/cluster.visualization
         done < ${curdir}/Analysis/2-MDS.coord

         #Plot  
         cd ${curdir}/Analysis/
         Rscript -e 'library(ggplot2); dd <- read.table("cluster.visualization", header = TRUE); dd$cluster <- factor(dd$cluster); p <- ggplot(data = dd, aes(x = M1,y = M2, color = cluster)) + geom_point(pch = 1, size = 3); ggsave("cluster_visualization.svg", width = 12, height = 12)'

         read -p $'\nCheck cluster_visualization.svg.\nAre you happy with the number of cluster?(y/n) ' vkyn
         if [ "${vkyn}" == "y" ] || [ "${vkyn}" == "Y" ]; then
         	((sc++))
         fi
      done
   fi
fi

((partindex++))
duration=$SECONDS
tstamp=$(date)
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${tstamp}
--------------------------------------------------------------------------------"
#*************************************************************

#5.Post analysis operation
#*************************************************************
echo "
Part ${partindex}. Post analysis operation"
SECONDS=0

#read -p "continue? " yn
yn="y"
if [ "${yn}" == "y" ]; then

if [ "${optik}" == "clustree" ]; then
   pknum=${vk}

   cd ${predir}
   #Create partition-sorted file
   rm -rf ${predir}/cl_sort
   mkdir ${predir}/cl_sort
   if [ "${pknum}" -gt 10 ]; then pknummd=$(echo "k${pknum}"); else pknummd=$(echo "k0${pknum}"); fi
   partSort ${predir}/clinfo/${pknummd}.cl
   mv ${predir}/clinfo/${pknummd}.cl_output ${predir}/cl_sort/

   #Execute the sequence splitting
   ni=0
   for file in ${predir}/*phy
   do       
     ((ni++))
     name=$(basename "${file}" | cut -d "." -f 1)
     mkdir ${predir}/${name}
     python ${amasdir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ${predir}/cl_sort/${pknummd}.cl_output -u phylip > /dev/null
     mv ${predir}/*cluster*phy ${predir}/${name}/
   done
   
   #Supergene tree construction
   for i in ${predir}/${name}/*phy
   do
     trename=$(basename "${i}" | cut -d "." -f 1 | cut -d "-" -f 1 | cut -d "_" -f 3)
     iqtree -s ${i} ${mf} -m ${kmodel} -ninit 200 -ntop 50 -nt 1 -seed ${iqseed} -keep-ident -quiet
     nw_reroot ${i}.treefile ZF > ${predir}/${name}/${trename}.tre
   done
   
   for t in ${predir}/${name}/*tre
   do
     sgtname=$(basename ${t} | cut -d "." -f 1)
     sgt=$(cat ${t})
     echo "${sgtname} ${sgt}" >> ${curdir}/gnTrees_collection.tre
     rm -rf ${t}
   done
fi

if [ -e ${predir}/logLkhforR.csv ]; then
#Find potential number of clusters
   cd ${predir}/
   Rscript <(echo "${find_potent_knum}") 2>&1 >/dev/null
   pknumall=$(cat potent_knum)
   pknum=$(sed -n "1p" potent_knum)
   read -p $"   
Potential number of clusters are
${pknumall}
Continue with ${pknum} clusters? (y/n): " pkyn
   #pkyn="y"
   if [ -z "${pkyn}" ] || [ "${pkyn}" == "Y" ] || [ "${pkyn}" == "y" ]; then 
      echo "
   *************************************************************************

               Msg:The analysis continues with ${pknum} clusters
     
   *************************************************************************"
   elif [ "${pkyn}" == "n" ] || [ "${pkyn}" == "N" ]; then
      read -p "Reset the number of clusters to: " pknumreset
      if [ -n "${pknumreset}" ] && [ "${pknumreset}" -eq "${pknumreset}" ] 2>/dev/null && [ "${pknumreset}" -le "${clfi}" ]; then
         pknum=${pknumreset}
         echo "
   *************************************************************************
     
              Msg:Number of clusters has been reset to ${pknum} 
                  and the analysis hence continues...
     
   *************************************************************************"
      else
         echo "
   *************************************************************************
     
              Msg:FAILED to reset...
                  the analysis continues with ${pknum} clusters
     
   *************************************************************************"
      fi
   else
      echo "
   *************************************************************************
     
              Msg:FAILED to reset...
                  the analysis continues with ${pknum} clusters
     
   *************************************************************************"
   fi

   if [ -n ${pknum} ]; then

      if [ "${pknum}" -gt 10 ]; then
         pknummd=$(echo "k${pknum}")
      else
         pknummd=$(echo "k0${pknum}")
      fi

      rm -rf ${curdir}/supergn_tree
      mkdir ${curdir}/supergn_tree
      
      mv ${predir}/AutoWZgn_concatenation/${pknummd}/*treefile ${curdir}/supergn_tree

      stepind=9
      #Extract the supergn_tree
      aftertreebuild ${curdir} supergn_tree ${stepind}
      mv ${curdir}/supergn_tree* ${curdir}/Trees/
     
   fi
fi

((partindex++))
duration=$SECONDS
tstamp=$(date)
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${tstamp}

--------------------------------------------------------------------------------"

fi
#*************************************************************

#6.Cluster Assignment
#*************************************************************
echo "
Part ${partindex}. Cluster assignment"
SECONDS=0

if [ -e ${curdir}/gnTrees_collection.tre ]; then
   #Distribute the supergene trees
   cd ${curdir}
   rm -rf supergn_trees
   mkdir supergn_trees   
   spi=1
   while [ "${spi}" -le "${pknum}" ]
   do
     name=$(cat gnTrees_collection.tre | sed -n "${spi}p" | cut -d " " -f 1)
     nwt=$(cat gnTrees_collection.tre | sed -n "${spi}p" | cut -d " " -f 2)     
     echo "${nwt}" > ./supergn_trees/${name}.tre
     ((spi++))
   done

   #Generate the tree collection file to calculate KCmetric distance
   rm -rf tree.KCmetric
   mkdir tree.KCmetric
   for sgt in ${curdir}/supergn_trees/*
   do
     sgtname=$(basename ${sgt} | cut -d "." -f 1)
     for bst in ${curdir}/bs-file_rooted/*
     do
       bstname=$(basename ${bst} | cut -d "." -f 1 | cut -d "-" -f 2)
       cat ${sgt} ${bst} > ${curdir}/tree.KCmetric/${bstname}-${sgtname}.tre
       echo "${curdir}/tree.KCmetric/${bstname}" >> ${curdir}/tree.KCmetric/kcdir-temp
     done
   done

   #Calculate the KCmetric distances between each supergene tree and the BS trees    
   cd ${curdir}/tree.KCmetric/
   KCmetric_temp=$(echo "${KCmetric}" | sed "s/args1/${lambda}/")
   Rscript <(echo "${KCmetric_temp}") 2>&1 >/dev/null
   rm -rf *tre
   for kdfile in *kcd
   do
     clsort=$(echo "${kdfile}" | cut -d "-" -f 1)
     if [ ! -d "${curdir}/tree.KCmetric/${clsort}" ]; then mkdir ${curdir}/tree.KCmetric/${clsort}; fi
     mv ${kdfile} ./${clsort}/
   done  

   #Assign clusters for each BS tree
   export -f MinFinder
   cat ${curdir}/tree.KCmetric/kcdir-temp | sort -n | uniq > ${curdir}/tree.KCmetric/kcdir
   rm -rf ${curdir}/tree.KCmetric/kcdir-temp
   cd ${curdir}/tree.KCmetric/
   parallel --no-notice -j ${thr} MinFinder :::: kcdir ::: ${curdir}/tree.KCmetric > /dev/null

   #Summary the information
   rm -rf ${curdir}/tree.KCmetric/Cluster_Assignment
   mkdir ${curdir}/tree.KCmetric/Cluster_Assignment

   fnum=$(ls $(sed -n "1p" ${curdir}/tree.KCmetric/kcdir) | wc -l)
   rm -rf ${curdir}/tree.KCmetric/kcdir
   
   if [ ! -e ${predir}/clinfo/order ]; then
   	  for i in ${curdir}/${alnf}/*phy
   	  do
   	  	name=$(basename "${i}" | cut -d "." -f 1)
        prefix=$(echo "${name}" | cut -d "-" -f 2 | cut -c 1-2)
        if [ "${prefix}" == "ns" ]; then
           echo "${name}" >> ${predir}/clinfo/orderpre.sex
        else
   	  	   echo "${name}" >> ${predir}/clinfo/orderpre.Auto
        fi
      done


      lnum=$(cat ${predir}/clinfo/orderpre.Auto | wc -l)
      l=1
      while [ "${l}" -le "${lnum}" ]
      do
        gene=$(cat ${predir}/clinfo/orderpre.Auto | sed -n "${l}p")
        line=$(cat ${predir}/clinfo/orderpre.Auto | sed -n "${l}p" | cut -d "-" -f 2)
        kindex=$(cat ${predir}/clinfo/${pknummd}.cl | grep "${line}" | cut -d " " -f 2)
        echo "${gene} ${kindex}" >> ${predir}/clinfo/ordersec.Auto
        ((l++))
      done
      rm -rf ${predir}/clinfo/orderpre.Auto
      
      cat ${predir}/clinfo/ordersec.Auto | sort -n -t ' ' -k 2 -k 1 | cut -d "-" -f 2 | cut -d " " -f 1 > ${predir}/clinfo/order.Auto
      rm -rf ${predir}/clinfo/ordersec.Auto

      if [ -e ${predir}/clinfo/orderpre.sex ]; then
   	     lnum=$(cat ${predir}/clinfo/orderpre.sex | wc -l)
   	     l=1
   	     while [ "${l}" -le "${lnum}" ]
   	     do
   	  	   gene=$(cat ${predir}/clinfo/orderpre.sex | sed -n "${l}p")
   	  	   line=$(cat ${predir}/clinfo/orderpre.sex | sed -n "${l}p" | cut -d "-" -f 2)
   	  	   kindex=$(cat ${predir}/clinfo/${pknummd}.cl | grep "${line}" | cut -d " " -f 2)
   	  	   echo "${gene} ${kindex}" >> ${predir}/clinfo/ordersec.sex
   	  	   ((l++))
   	     done
   	     rm -rf ${predir}/clinfo/orderpre.sex
      
         cat ${predir}/clinfo/ordersec.sex | sort -t ' ' -k 1 -k 2 | cut -d "-" -f 2 | cut -d " " -f 1 > ${predir}/clinfo/order.sex
         rm -rf ${predir}/clinfo/ordersec.sex

         cat ${predir}/clinfo/order.Auto ${predir}/clinfo/order.sex > ${predir}/clinfo/order
         rm -rf ${predir}/clinfo/order.Auto ${predir}/clinfo/order.sex
      else
         mv ${predir}/clinfo/order.Auto ${predir}/clinfo/order
         rm -rf ${predir}/clinfo/order.Auto
      fi
   fi


   odline=$(cat ${predir}/clinfo/order | wc -l)
   odi=1
   while [ "${odi}" -le "${odline}" ]
   do
     gnord=$(cat ${predir}/clinfo/order | sed -n "${odi}p")
     kcdasn=$(echo "${curdir}/tree.KCmetric/${gnord}-cluster.assign")   
     clori=${odi}
     if [ "${clori}" -lt 10 ]; then clori=$(echo "0${clori}"); fi
     name=$(echo -e "${clori}\t${gnord}")

     k=1
     while [ "${k}" -le "${fnum}" ]
     do
       kcount[$k]=0
       while read line
       do
         if [ "${line}" -eq "${k}" ]; then
            kcount[$k]=$((kcount[$k] + 1))
         fi
       done < ${kcdasn}

       if [ -z "${kcount[$k]}" ]; then kcount[$k]=0; fi       

       ((k++))
     done

     #Produce short-data table output
     l=2
     ksum=${kcount[1]}
     while [ "${l}" -le "${fnum}" ]
     do
       ksum=$(echo -e "${ksum}\t${kcount[$l]}")
       ((l++))
     done
     echo -e "${name}\t${ksum}" >> ${curdir}/tree.KCmetric/Cluster_Assignment_shortdata.summary        

     #Produce long-data table output
     l=1
     while [ "${l}" -le "${fnum}" ]
     do
       echo -e "${name}\tcluster${l}\t${kcount[$l]}" >> ${curdir}/tree.KCmetric/Cluster_Assignment.summary
       ((l++))
     done
     
     #File sorting
     mv ${curdir}/tree.KCmetric/${gnord}-cluster.assign ${curdir}/tree.KCmetric/Cluster_Assignment

     ((odi++))
   done

   #File sorting   
   rm -rf ${curdir}/tree.KCmetric/KCmetric_distance
   mkdir ${curdir}/tree.KCmetric/KCmetric_distance
   for alnfile in ${curdir}/${alnf}/*phy
   do
     alnname=$(basename ${alnfile} | cut -d "." -f 1 | cut -d "-" -f 2)
     mv ${curdir}/tree.KCmetric/${alnname} ${curdir}/tree.KCmetric/KCmetric_distance/
   done
   
   #Plot graphs for cluster assignment
   cd ${curdir}/tree.KCmetric/
   clusterAssign_temp=$(echo "${clusterAssign}" | sed "s/args1/${pknum}/")
   Rscript <(echo "${clusterAssign_temp}")  2>&1 >/dev/null
   
fi

duration=$SECONDS
echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.

--------------------------------------------------------------------------------"
#*************************************************************
rm -rf ${curdir}/supergn_trees 
mv ${curdir}/gnTrees_collection.tre ${curdir}/supergnTrees_collection.tre
mv ${curdir}/bs-file_rooted ${curdir}/supergnTrees_collection.tre ${curdir}/Trees/

tstamp=$(date)
echo "
Congratulations!!
Analysis finished at ${tstamp}"
fi
