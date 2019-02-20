#!/bin/bash

######################
# Function 1. Manual
######################
usage() {
  echo "Usage: [-a] alignment folder
       [-b] trees folder
       [-c] use codon model
       [-d] distance metirc: "gd\(Geodesic\)", "kc\(Kendall-Colijn\)", "pd\(Probabilistic distance\)"
       [-h] show usage
       [-L] a lambda between 0 and 1; NA in Geodesic and Probabilistic distance
       [-m] cluster method: "ward\(Default\)", "kmeans", "MDSK", "pam", "MDSP", "NMDSK"
       [-k] infer optimal number of clusters: "gori\(Default\)", "clustree"
       [-o] specify outgroup
       [-p] "f" to use fast method; "s" to use slow method; only valid if [-c]
       [-t] number of threads
       [-s] seed number for all tree construction"
}

#########################################################
# Function 2. Extract and root the trees after building
#########################################################
AfterTreesConstruction() {
  rm -rf ${1}/${2}_rooted
  mkdir ${1}/${2}_rooted
  for file in ${1}/${2}/*treefile; do
    name_raw=$(echo "${file}")
    name=$(basename ${name_raw})
    nw_reroot ${file} ${3} >${1}/${2}_rooted/${name}
  done
  if [[ "${4}" -gt 1 ]]; then
    for file in ${1}/${2}_rooted/*treefile; do
      name_raw=$(echo "${file}")
      if [[ "${4}" -gt 5 ]]; then
        name=$(basename ${name_raw} | cut -d "." -f 1 | cut -d "_" -f 3 | cut -d "-" -f 1)
      else
        name=$(basename ${name_raw} | cut -d "." -f 1 | cut -d "-" -f 2)
      fi
      seq=$(sed -n "1p" ${file})
      echo "$name $seq" >>${1}/gnTrees_collection.tre
    done
  fi
}

####################################################################
# Function 3. Make permutation datasets from the output of fseqboot
####################################################################
MakeDatasets() {
  echo -e "The seed for the generation of the 100-replicate permutation data is: ${3}." >impMsg.txt
  i=1
  j=$((i + 100))
  file_number=1
  percent=5
  rm -rf ${1}/permutation_data
  mkdir ${1}/permutation_data
  until [[ "${i}" -eq "${j}" ]]; do
    k=$((i + 1))
    start_number=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${i}p")
    end=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${k}p")
    end_number=$((end - 1))
    cat ${1}/outfile | sed -n "${start_number},${end_number}p" >${1}/permutation_data/d${file_number}.phy
    i=$((i + 1))
    file_number=$((file_number + 1))
  done
}

###################################################
# Function 4. Generate file containing partitions
###################################################
GeneratePartitionFile() {
  seq_pos1=1
  seq_tlen=0
  for file in ${1}/${2}/*phy; do
    name=$(basename "${file}" | cut -d "." -f 1 | cut -d "-" -f 2)
    seq=$(cat ${file} | cut -d " " -f 3)
    seq_tlen=$(echo "$((seq_tlen + seq))")
    echo "${name} = $seq_pos1-$seq_tlen;" >>${1}/Analysis/${3}/partition
    seq_pos1=$(echo "$((seq_tlen + 1))")
  done
}

##############################################################################
# Function 5. Sorting out the partitions according to the clustering results
##############################################################################
PartitionSorting() {
  output_name=$(echo "${1}_output")
  line_number=$(cat partition | wc -l)
  i=1
  while [[ "${i}" -le "${line_number}" ]]; do
    gene_name[$i]=$(cat partition | sed -n "${i}p" | cut -d " " -f 1)
    partition[$i]=$(cat partition | sed -n "${i}p" | cut -d " " -f 3 | cut -d ";" -f 1)
    ((i++))
  done
  i=1
  while [[ "${i}" -le "${line_number}" ]]; do
    cluster[$i]=$(cat ${1} | grep "${gene_name[$i]}" | cut -d " " -f 2)
    ((i++))
  done
  rm -f ${output_name}
  i=1
  k=1
  while [[ "${i}" -le "${line_number}" ]]; do
    if [[ "${i}" -eq "${line_number}" ]]; then
      j=${i}
    else
      j=$((i + 1))
    fi
    if [[ "${cluster[$i]}" != "NA" ]]; then
      partition_content=$(echo "${partition[$i]}")
      while [[ "${j}" -le "${line_number}" ]]; do
        if [[ "${cluster[$j]}" != "NA" ]]; then
          if [[ "${cluster[$j]}" -eq "${cluster[$i]}" ]]; then
            if [[ "${i}" -lt "${line_number}" ]]; then
              partition_content=$(echo "${partition_content},${partition[$j]}")
              cluster[$j]="NA"
            fi
          fi
        fi
        ((j++))
      done
      if [[ "${k}" -lt 10 ]]; then
        echo "cluster0${k} = ${partition_content}" >>${output_name}
      else
        echo "cluster${k} = ${partition_content}" >>${output_name}
      fi
      ((k++))
    fi
    ((i++))
  done
}

##############################################################
# Function 6. Extract and compute log likelihood of datasets
##############################################################
loglkhExtraction() {
  cluster_max_file_number=$(ls ${2}/*iqtree | wc -l)
  k=1
  while [[ "${k}" -le "${cluster_max_file_number}" ]]; do
    if [[ "${k}" -lt 10 ]]; then
      k_mod=$(echo "0${k}")
    else
      k_mod=${k}
    fi
    loglkh=$(cat ${2}/*cluster${k_mod}*iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
    loglkh_mod=$(echo "${loglkh} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
    loglkh_mod_sum=$((loglkh_mod_sum + loglkh_mod))
    ((k++))
  done
  seq_name=$(echo "${2}" | rev | cut -d "/" -f 2 | rev)
  echo -e "${seq_name} k=${cluster_max_file_number} -loglikelihood*10^4:\t${loglkh_mod_sum}" >>${1}/${seq_name}_logLhcollection
}

#######################################################
# Function 7. A Fast method to permute based on codon
#######################################################
CodonJackniffPermutation() {
  codon_permute_dir=${1}
  len=${2}
  d=${3}
  codon_sum=$((len / 3))
  f=$(((codon_sum / 32767) + 1))
  i=1
  until [[ "${i}" -gt 1 ]]; do
    r=$(((${RANDOM} * f) / 2 + (codon_sum / 10)))
    if [[ "${r}" -lt "$((codon_sum / 2))" ]]; then
      iteration=${r}
      ((i++))
    fi
  done
  cp -p ${codon_permute_dir}/aln.ori ${codon_permute_dir}/aln.ori-${d}
  i=1
  while [[ "${i}" -le "${iteration}" ]]; do
    ri=1
    until [[ "${ri}" -gt 1 ]]; do
      r=$(((${RANDOM} * f) + 1))
      if [[ "${r}" -lt "${codon_sum}" ]]; then
        p=${r}
        ((ri++))
      fi
    done
    permute_position=$((p * 3))
    cat ${codon_permute_dir}/aln.ori-${d} | cut -c 1-${permute_position} >${codon_permute_dir}/${d}p1.txt
    cat ${codon_permute_dir}/aln.ori-${d} | cut -c $((permute_position + 1))- >${codon_permute_dir}/${d}p2.txt
    paste -d "" ${codon_permute_dir}/${d}p2.txt ${codon_permute_dir}/${d}p1.txt >${codon_permute_dir}/aln.ori-${d}0
    mv ${codon_permute_dir}/aln.ori-${d}0 ${codon_permute_dir}/aln.ori-${d}
    rm -f ${codon_permute_dir}/aln.ori-${d}0 ${codon_permute_dir}/${d}p1.txt ${codon_permute_dir}/${d}p2.txt
    ((i++))
  done
  paste -d "" ${codon_permute_dir}/taxa ${codon_permute_dir}/aln.ori-${d} >${codon_permute_dir}/d${d}.part
  cat ${codon_permute_dir}/info ${codon_permute_dir}/d${d}.part >${codon_permute_dir}/d${d}.phy
  rm -f ${codon_permute_dir}/aln.ori-${d} ${codon_permute_dir}/d${d}.part
}

######################################
# Function 8. Probabilistic distance
######################################
ProbabilisticDistanceCalculator() {
  dir=${1}
  distance_type=$(echo "-${2}")
  exact_distance=${3}
  # Approximation type
  if [[ "${2}" == "h" ]]; then
    approx_type="-ecv"
  else
    approx_type="-e"
  fi
  # Tree file modification
  cat ${dir}/Analysis/gnTrees_collection.tre | cut -d " " -f 1 >tree_name.list
  cat ${dir}/Analysis/gnTrees_collection.tre | cut -d " " -f 2 >gnTrees_collection_mod.tre
  # Generate a file supplement the self-comparison distance
  tree_number=$(cat ${dir}/Analysis/gnTrees_collection.tre | wc -l)
  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    echo -e "\nTree ${i}, Tree ${i}, Alignment length NA, Distance 0\n" >>${dir}/Analysis/sc_dist.prob
    ((i++))
  done
  # Calculate Probabilistic Distance
  redo_calculation=0
  while [[ "${redo_calculation}" -lt 1 ]]; do
    if [[ "${distance_type}" == "-h" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Hellinger distance." >>${work_dir}/parameter_input
    elif [[ "${distance_type}" == "-kl" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Kullback-Leibler divergence." >>${work_dir}/parameter_input
    elif [[ "${distance_type}" == "-js" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Jensen-Shannon distance." >>${work_dir}/parameter_input
    fi
    # Calculation begins
    if [[ "${exact_distance}" == "ex" ]]; then
      if [[ "${redo_calculation}" -eq 0 ]]; then
        echo "   Probabilistic distance is calcalated based on trees only." >>${work_dir}/parameter_input
      fi
      java -classpath "${probdist_dir}/*" simplemetrics.ProbabilityDistanceApp -ex ${distance_type} ${dir}/Analysis/gnTrees_collection_mod.tre >${dir}/Analysis/tree_dist_raw.prob
    else
      java -classpath "${probdist_dir}/*" simplemetrics.ProbabilityDistanceApp -m ${dir}/Analysis/modelParamFile.txt ${distance_type} ${approx_type} -re 0.05 -zb 1.645 -aug ${dir}/Analysis/gnTrees_collection_mod.tre >${dir}/Analysis/tree_dist_raw.prob
    fi
    na_exist=$(cat ${dir}/Analysis/tree_dist_raw.prob | grep "NaN")
    if [[ -z "${na_exist}" ]]; then
      ((redo_calculation++))
    else
      echo -e "\nNaN occurred! Probabilistic distance would be re-calculated."
    fi
  done
  # File modification
  cat tree_dist_raw.prob sc_dist.prob >tree_dist.prob
  rm -f tree_dist_raw.prob sc_dist.prob gnTrees_collection_mod.tre
  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    if [[ "${i}" -lt 10 ]]; then
      i_mod=$(echo "Tree0${i}")
    else
      i_mod=$(echo "Tree${i}")
    fi
    sed -i "s/Tree ${i},/${i_mod}/g" tree_dist.prob
    ((i++))
  done
  Rscript <(echo "${ProbDistMatrix}") 2>&1 >/dev/null
  #rm -f tree_dist.prob
  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    if [[ "${i}" -lt 10 ]]; then
      i_mod=$(echo "Tree0${i}")
    else
      i_mod=$(echo "Tree${i}")
    fi
    tree_name=$(cat ${dir}/Analysis/tree_name.list | sed -n "${i}p")
    sed -i "s/${i_mod}/${tree_name}/g" ${dir}/Analysis/dist_matrix.txt
    ((i++))
  done
  rm -f tree_name.list
  #Clustering
  Rscript <(echo "${ClusteringProbdist}") 2>&1 >/dev/null
  # File sorting
  #rm -f ${dir}/Analysis/dist_matrix.txt
}

#############################
# Function 9. Phylip format
#############################
PhylipFormat() {
  dir=${1}
  line_number=$(cat ${dir}/aln_concatenation.phy | wc -l)
  line_index=2
  while [[ "${line_index}" -le "${line_number}" ]]; do
    col=1
    while [[ "${col}" -le 25 ]]; do
      a=$(cat ${dir}/aln_concatenation.phy | sed -n "${line_index}p" | cut -c ${col})
      if [[ "${a}" == " " ]]; then
        col_next=$((col++))
        b=$(cat ${dir}/aln_concatenation.phy | sed -n "${line_index}p" | cut -c ${col_next})
        if [[ "${b}" == " " ]]; then
          ai=${col}
          col=$((col + 25))
        fi
      fi
      ((col++))
    done
    space_number=$(cat ${dir}/aln_concatenation.phy | sed -n "${line_index}p" | cut -c ${ai}-25 | grep -o ' ' | wc -l)
    space_to_fill=$((10 - space_number - ai + 1))
    s=1
    space=" "
    while [[ "${s}" -lt "${space_to_fill}" ]]; do
      space=$(echo "${space} ")
      ((s++))
    done
    taxa_name_end_col=$((ai - 1))
    taxa_name=$(cat ${dir}/aln_concatenation.phy | sed -n "${line_index}p" | cut -c 1-${taxa_name_end_col})
    sed -i "s/${taxa_name}/${taxa_name}${space}/" ${dir}/aln_concatenation.phy
    ((line_index++))
  done

}

##########################
# R code 1. Get clusters
##########################
Clustering='if (!require(treespace)) install.packages("treespace")
if (!require(heatmap3)) install.packages("heatmap3")
if (!require(cluster)) install.packages("cluster")
if (!require(smacof)) install.packages("smacof")
if (!require(distory)) install.packages("distory")
if (!require(phangorn)) install.packages("phangorn")

trees <- read.tree("gnTrees_collection.tre")
cluster_number <- length(trees)

if ("args1" == "KC_distance") {
	L <- args2
    KC_res <- treespace(trees, method = "treeVec", nf = 2, lambda = L, return.tree.vectors = TRUE)
    distance <- KC_res$D
} else if ("args1" == "Geodesic_distance") {
	geod <- dist.multiPhylo(trees)
    wrfd <- wRF.dist(trees)

    wrfd <- as.matrix(wrfd)
    treename<- rownames(wrfd)
    geod <- as.matrix(geod)
    rownames(geod) <- treename
    colnames(geod) <- treename
    distance <- as.dist(geod)
}

if ("args3" == "ward") {
   ##################################################
   # Hierachichal clustering with MDS visualization
   ##################################################
   mdsres2D <- smacofSym(distance, ndim = 2, type = "ratio")
   treMDS2D <- as.data.frame(mdsres2D$conf)
   write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
   write.table(mdsres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
   mdsres3D <- smacofSym(distance, ndim = 3, type = "ratio")
   treMDS3D <- as.data.frame(mdsres3D$conf)
   write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
   write.table(mdsres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)
   distance.hclust <- hclust(distance, method = "ward.D2")
   svg("dendrogram.svg", width=10, height=10)
   plot(distance.hclust, cex = 0.6)
   dev.off()
   distance_matrix <- as.matrix(distance)
   svg("heatmap.svg", width=12, height=12)
   heatmap3(distance_matrix, Rowv = as.dendrogram(distance.hclust), symm = TRUE)
   dev.off()
   dend_label_raw <- as.data.frame(distance.hclust$labels)
   dend_label <- dend_label_raw$`distance.hclust$labels`[order.dendrogram(as.dendrogram(distance.hclust))]
   dend_label <- as.data.frame(dend_label)
   write.table(dend_label,file = "dend.order", quote = FALSE, col.names = FALSE, row.names = FALSE)
   if (cluster_number > 15){cluster_number<-15}
   for(j in 1:cluster_number){
      cl <- cutree(distance.hclust, k = j)
       if (j < 10){
          j=paste("0",j,sep= "")
          nj=paste("k",j,".cl",sep = "")
        } else {
          nj=paste("k",j,".cl",sep = "")
        }
      write.table(cl, nj, quote = FALSE, col.names = FALSE)
    }
}
###########################
# K-means on tree vectors
###########################
if (("args1" == "KC_distance") & ("args3" == "kmeans")) {
    treMDS <- KC_res$pco$li
    write.table(treMDS, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    treVec <- KC_res$vectors
    if (cluster_number > 15){
   	  cluster_number <- 15
    } else {
      cluster_number <- cluster_number - 1
    }
    for (j in 1: cluster_number){
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
##################################################################
# Multi-Dimensional scaling down to 2D and 3D, K-means on 3D MDS
##################################################################
if ("args3" == "MDSK") {
    mdsres2D <- smacofSym(distance, ndim = 2, type = "ratio")
    treMDS2D <- as.data.frame(mdsres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdsres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    mdsres3D <- smacofSym(distance, ndim = 3, type = "ratio")
    treMDS3D <- as.data.frame(mdsres3D$conf)
    write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdsres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)
    if (cluster_number > 15){
       cluster_number <- 15
    } else {
       cluster_number <- cluster_number - 1
    }
    for (j in 1: cluster_number){
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
##############################################
# PAM on the distance/dissimilarities matrix
##############################################
if ("args3" == "pam") {
    if (cluster_number > 15){
   	  cluster_number <- 15
    } else {
      cluster_number <- cluster_number - 1
    }
    for (j in 1: cluster_number){
        kc.pam <- pam(distance, diss = TRUE, k = j, cluster.only = TRUE)
        if (j < 10){
          j=paste("0",j,sep= "")
          nj=paste("k",j,".cl",sep = "")
        } else {
          nj=paste("k",j,".cl",sep = "")
        }
        write.table(kc.pam, nj, quote = FALSE, col.names = FALSE)
    }
}
##############################################################
# Multi-Dimensional scaling down to 2D and 3D, PAM on 3D MDS
##############################################################
if ("args3" == "MDSP") {
    mdsres2D <- smacofSym(distance, ndim = 2, type = "ratio")
    treMDS2D <- as.data.frame(mdsres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdsres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    mdsres3D <- smacofSym(distance, ndim = 3, type = "ratio")
    treMDS3D <- as.data.frame(mdsres3D$conf)
    write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdsres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)
    if (cluster_number > 15){
       cluster_number <- 15
    } else {
       cluster_number <- cluster_number - 1
    }
    for (j in 1: cluster_number){
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
###############################################################
# Non-metric Multi-Dimensional scaling down to 2D and K-means
###############################################################
if ("args3" == "NMDSK") {
    mdsres2D <- smacofSym(distance, ndim = 2, type = "ordinal")
    treMDS2D <- as.data.frame(mdsres2D$conf)
    write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
    write.table(mdsres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
    if (cluster_number > 15){
       cluster_number <- 15
    } else {
       cluster_number <- cluster_number - 1
    }
    for (j in 1: cluster_number){
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

################################################################
# R code 2. Format table for visualization of permutation test
################################################################
loglkhImprovementMake='out.file <- ""
file.names <- dir(pattern = "logLhcollection")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i])
  out.file <- cbind(out.file, file)
}
write.csv(out.file, file = "logLkhcollection.csv")
da <- read.csv("logLkhcollection.csv")
actual_data <- read.table("aln_concatenation_logLkhcollection", header = FALSE)
da <- cbind(actual_data,da[,3:102])
permutation_data_number <- seq(1:100)
names(da) <- c("Actual.data", permutation_data_number)
rownum <- nrow(da)
for(i in 1:rownum){
	da[i,] <- 0 - (da[i,] / 10000)
}
permutation_data_header <- rep(" ", rownum)
kindex <- seq(1: rownum)
out.file <- cbind(kindex,da$Actual.data,permutation_data_header,da[,2:101])
names(out.file) <- c("k","Actual data","Permutation data",permutation_data_number)
write.csv(out.file, file = "logLkhcollection.csv", row.names = FALSE)
for(k in 1:(rownum - 1)){
  da[k,] <- da[k+1,] - da[k,]
}
da <- rbind(da[1:(rownum - 1),])
dkindex <- seq(1:(rownum - 1))
permutation_data_mean <- apply(cbind(da[,2:101]), 1, mean)
permutation_data_sd <- apply(cbind(da[,2:101]), 1, sd)
daforR <- cbind(dkindex,da$Actual.data,permutation_data_mean,permutation_data_sd,da[,2:101])
names(daforR) <- c("dk","Actual data","Average","SD",permutation_data_number)
write.csv(daforR,file = "logLkhforR.csv", row.names = FALSE)'

#############################################
# R code 3. Plot graph for permutation test
#############################################
loglkhImprovementPlot='if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
loglkh_data <- read.csv("logLkhforR.csv")
a <- loglkh_data[,2:3]
write.csv(a,file = "a.csv")
a <- read.csv("a.csv")
a$X <- c("k01","k02","k03","k04","k05","k06","k07","k08","k09","k10","k11","k12","k13","k14")
a <- as.data.frame(a)
names(a) <- c("cluster_num","Actual.data","Average")
loglkh_data <- t(loglkh_data)
b <- loglkh_data[5:104,]
write.csv(b, file = "b.csv")
b <- read.csv("b.csv")
names(b) <- c("X","k01","k02","k03","k04","k05","k06","k07","k08","k09","k10","k11","k12","k13","k14")
num <- seq(1:100)
b$X <- num
bb <- b[,2:14]
write.csv(bb, file = "b.csv")
a <- gather(a,group,dk,Actual.data:Average)
b <- gather(b,cluster_num,dk,k01:k14)
write.csv(a,file = "amd.csv")
write.csv(b,file = "bmd.csv")
a <- read.csv("amd.csv")
a <- as.data.frame(a)
b <- read.csv("bmd.csv")
b <- as.data.frame(b)
xlabel <- c("1→2", "2→3", "3→4", "4→5", "5→6", "6→7", "7→8", "8→9", "9→10", "10→11", "11→12", "12→13", "13→14", "14→15")
bp <- ggplot(b,aes(cluster_num,dk)) + geom_boxplot(data=b,aes(cluster_num,dk), color="grey") + geom_point(data=a, aes(x=cluster_num,y=dk, group=group, color = group)) + geom_line(data=a, aes(x=cluster_num,y=dk, group = group, color = group)) + scale_x_discrete(labels = xlabel) + labs(x = "Number of clusters, k", y = "log Likelihood improvement") + theme(axis.text.x = element_text(size = 18, vjust = 0.5, angle = 45, family = "Calibri"), axis.title.x = element_text(size = 20, family = "Calibri"), axis.text.y = element_text(size = 18, family = "Calibri"), axis.title.y = element_text(size = 20, family = "Calibri"), legend.title = element_text(size = 20, family = "Calibri"), legend.text = element_text(size = 20, family = "Calibri")) + scale_color_manual(name = "Legend", values = c("firebrick2", "royalblue2"), labels = c("Actual data", "Permutation data mean"))
pdf(NULL)
ggsave("permutation_test.svg",bp, width = 14, height = 10)'

####################################################
# R code 4. A way to detect optimal cluster number
####################################################
InferPotentialk='loglkh_data <- read.csv("logLkhforR.csv", header = TRUE)
loglkh_data <- abs(loglkh_data[,2]-loglkh_data[,3])
loglkh_min <- min(loglkh_data, na.rm = TRUE)
min_index <- which(loglkh_data == loglkh_min)
loglkh_percent <- seq(1:14)
for (i in 1:13) {
  loglkh_percent[i+1] <- loglkh_data[i+1]/loglkh_data[1]
}
res <- array("NA", 14)
j <- 1
threshd <- 0.5
for (i in 1:12){
  if (i >= 2){
    if (loglkh_percent[i] <= threshd){
      a <- loglkh_data[i-1] - loglkh_data[i]
      b <- loglkh_data[i] - loglkh_data[i+1]
      if (b <= a) {
        if (loglkh_percent[i+1] <= threshd){
        res[j] <- i
        j <- j+1
        }
      }
    }
  }
  else {
    if (loglkh_data[i] == loglkh_min){
      res[j] <- i
      j <- j+1
    }
  }
}
res <- res[res != "NA"]
res <- as.data.frame(res)
write.table(res, file = "potent_knum", quote = FALSE, col.names = FALSE, row.names = FALSE)'

######################################################################################
# R code 5. Compute square of difference between alternative model and optimal model
######################################################################################
dBICCalculator='if (!require(foreach)) install.packages("foreach")
if (!require(doParallel)) install.packages("doParallel")
number_cores <- detectCores()
registerDoParallel(number_cores)
files <- list.files(pattern = "*model")
foreach(i=files) %dopar% {
    name <- strsplit(i, ".", fixed = TRUE)[[1]][1]
    name <- paste(name, ".txt", sep = "")
  dd <- read.table(i)
  optimal_model <- dd[1,2]
  for (k in 1:nrow(dd)){
    dd[k,2] <- round(((dd[k,2] - optimal_model) ^ 2), 3)
  }
  write.table(dd[order(dd$V1),], file = name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}'

########################################################
# R code 6. Find the single model for permutation test
########################################################
InferSingleModel='if (!require(dplyr)) install.packages("dplyr")
dd <- read.table("model.summary")
dd$sum_sq <- apply(dd[,2:ncol(dd)], 1, sum)
dd <- arrange(dd, sum_sq)
write.table(dd, "bestmodel.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
'

###########################
# R code 7. Permute array
###########################
ArrayPermutation='array <- read.table("pos.array")
seed_number <- seednumber
set.seed(seed_number)
for (i in 1:5){
    seed_number <- seed_number + 1
    set.seed(seed_number)
	p <- sample(array)
	name <- paste(i,".shuffle", sep="")
	write.table(p, file = name, quote = FALSE, row.names = FALSE, col.names = FALSE)
	}
'

##################################
# R code 8. Make ProbDist matrix
##################################
ProbDistMatrix='if (!require(reshape2)) install.packages("reshape2")
dd <- read.table("tree_dist.prob")
dd <- data.frame(dd$V1,dd$V2,dd$V7)
names(dd) <- c("a","b","value")
dm <- as.dist(acast(dd,b~a, value.var = "value"))
dm <- as.matrix(dm)
write.table(dm,"dist_matrix.txt", quote = F)
'

###############################################################
# R code 9. Hierarchical clustering on Probabilistic distance
###############################################################
ClusteringProbdist='if (!require(ape)) install.packages("ape")
if (!require(heatmap3)) install.packages("heatmap3")
if (!require(cluster)) install.packages("cluster")
if (!require(smacof)) install.packages("smacof")

dm <- read.table("dist_matrix.txt", row.names = 1)
dm <- as.dist(dm)
###########################
# Hierachichal Clustering
###########################
trees <- read.tree("gnTrees_collection.tre")
cluster_number <- length(trees)
mdsres2D <- smacofSym(dm, ndim = 2, type = "ratio")
treMDS2D <- as.data.frame(mdsres2D$conf)
write.table(treMDS2D, file = "2-MDS.coord", col.names = FALSE, quote = FALSE)
write.table(mdsres2D$stress, file = "2-MDS.stress", col.names = FALSE, quote = FALSE)
mdsres3D <- smacofSym(dm, ndim = 3, type = "ratio")
treMDS3D <- as.data.frame(mdsres3D$conf)
write.table(treMDS3D, file = "3-MDS.coord", col.names = FALSE, quote = FALSE)
write.table(mdsres3D$stress, file = "3-MDS.stress", col.names = FALSE, quote = FALSE)
dm.hclust <- hclust(dm, method = "ward.D2")
svg("dendrogram.svg", width=10, height=10)
plot(dm.hclust, cex = 0.6)
dev.off()
dm <- as.matrix(dm)
svg("heatmap.svg", width=12, height=12)
heatmap3(dm, Rowv = as.dendrogram(dm.hclust), symm = TRUE)
dev.off()
dend_label_raw <- as.data.frame(dm.hclust$labels)
dend_label <- dend_label_raw$`dm.hclust$labels`[order.dendrogram(as.dendrogram(dm.hclust))]
dend_label <- as.data.frame(dend_label)
write.table(dend_label,file = "dend.order", quote = FALSE, col.names = FALSE, row.names = FALSE)
if (cluster_number > 15){cluster_number<-15}
for(j in 1:cluster_number){
    cl <- cutree(dm.hclust, k = j)
    if (j < 10){
        j=paste("0",j,sep= "")
        nj=paste("k",j,".cl",sep = "")
    } else {
        nj=paste("k",j,".cl",sep = "")
    }
    write.table(cl, nj, quote = FALSE, col.names = FALSE)
}
'

###################################################
# Pre-step: Load the flags and check the settings
###################################################
while getopts "ha:b:cd:u:L:m:k:o:p:s:t:" opt; do
  case $opt in
  h) usage && exit ;;
  a) aln_folder=$(echo "${OPTARG}" | cut -d "/" -f 1) ;;
  b) tree_folder=$(echo "${OPTARG}" | cut -d "/" -f 1) ;;
  c) codon_model="Yes" ;;
  d) distance_metric=${OPTARG} ;;
  L) lambda=${OPTARG} ;;
  m) cluster_method=${OPTARG} ;;
  k) infer_optimal_k=${OPTARG} ;;
  o) out_group=${OPTARG} ;;
  p) permutation_method=${OPTARG} ;;
  s) iqtree_seed_number=${OPTARG} ;;
  t) thread=${OPTARG} ;;
  esac
done

# No argument returns usage manual and exit
if [[ $OPTIND -eq 1 ]]; then
  usage
  exit
fi

# Get current directory
work_dir=$(pwd)

# Check if alignment folder is specified correctly
if [[ -z "${aln_folder}" ]] || [[ ! -d "${work_dir}/${aln_folder}" ]]; then
  echo -e "\nMsg: Some files are missing.\n"
  usage
  exit
fi

# Check if trees folder is specified
if [[ -z "${tree_folder}" ]]; then
  tree_folder="NA"
fi

# Check if outgroup is specified
if [[ -z "${out_group}" ]]; then
  file=$(ls ${work_dir}/${aln_folder} | sed -n "1p")
  out_group=$(cat ${work_dir}/${aln_folder}/${file} | sed -n "2p" | cut -d " " -f 1)
fi

# Check if the distance metric is specified
if [[ -z "${distance_metric}" ]]; then
  distance_metric="kc"
fi

# Check if codon model is used
if [[ "${codon_model}" == "Yes" ]]; then
  model_number=92
  codon_model_flag="-st CODON"
  evo_model_display="Codon substitution models"
  if [[ "${distance_metric}" == "pd" ]]; then
    echo -e "\nMsg: Probabilistic distance is NOT applicable when codon model is chosen."
    exit
  fi
else
  model_number=134
  codon_model_flag=""
  evo_model_display="Nucleotide substitution models"
fi

# Check if probabilistic distance is used
if [[ "${distance_metric}" == "pd" ]] || [[ "${distance_metric}" == "pd_ex" ]]; then
  distance_metric_display="Probabilistic distance"
  lambda="NA"
elif [[ "${distance_metric}" == "gd" ]]; then
  distance_metric_display="Geodesic distance"
  lambda="NA"
elif [[ "${distance_metric}" == "kc" ]]; then
  distance_metric_display="Kendall-Colijn metric"
fi

# Check the lambda setting
if [[ "${distance_metric}" == "kc" ]]; then
  if [[ -z "${lambda}" ]]; then
    echo -e "\nERROR: Missing lambda value!"
    exit
  fi
fi

# Check the clustering method setting
# Default: Hierachichal clustering with Ward.D2 linkage
if [[ -z "${cluster_method}" ]]; then
  cluster_method="ward"
fi

# Cluster method info
case ${cluster_method} in
ward) cluster_method_display="Hierachichal clustering" ;;
kmeans) cluster_method_display="k-means clustering" ;;
MDSK) cluster_method_display="MDS -> k-means clustering" ;;
pam) cluster_method_display="k-medoids (PAM) clustering" ;;
MDSP) cluster_method_display="MDS -> k-medoids (PAM) clustering" ;;
NMDSK) cluster_method_display="NMDS" ;;
esac

# Check method for inferring optimal cluster number
if [[ -z "${infer_optimal_k}" ]]; then
  case ${cluster_method} in
  ward) infer_optimal_k="gori" ;;
  kmeans) infer_optimal_k="clustree" ;;
  MDSK) infer_optimal_k="clustree" ;;
  pam) infer_optimal_k="gori" ;;
  MDSP) infer_optimal_k="clustree" ;;
  NMDSK) infer_optimal_k="clustree" ;;
  esac
fi

# If number of thread isn't specified
# detect the number of processors
if [[ -z "${thread}" ]]; then
  thread=$(nproc --all)
fi

# If seed isn't specified for IQ-Tree
# generate a random seed
if [[ -z "${iqtree_seed_number}" ]]; then
  iqtree_seed_number=${RANDOM}
fi

# Resize the terminal window
resize -s 36 80 >/dev/null

# Show the requirements and settings
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
read -p " Files checklist: 1.Alignment folder

    Dependencies: 1.IQ-Tree               4.fseqboot
                  2.Newick Utilities      5.probdist
                  3.AMAS                  6.Rstudio

      R packages: 1.treespace              6.smacof
                  2.heatmap3               7.dplyr
                  3.ggplot2                8.foreach
                  4.reshape2               9.doParallel
                  5.tidyr                 10.clustree
   *****************************************************************
      SETTINGS:
                 Alignment folder: ${aln_folder}
                     Trees folder: ${tree_folder}
                         Outgroup: ${out_group}
                  Evolution model: ${evo_model_display}
                  Distance metric: ${distance_metric_display}
                           Lambda: ${lambda}
                   Cluster method: ${cluster_method_display}
                  Infer optimal k: ${infer_optimal_k}
                          Threads: ${thread}
                             Seed: ${iqtree_seed_number}
   *****************************************************************

Continue?(Y/n): " answer

####################
# Programme starts
####################
if [[ "${answer}" == "Y" ]] || [[ "${answer}" == "y" ]]; then

  # Generate output containing settings
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
                /_/
                ' >${work_dir}/parameter_input
  echo "
SETTINGS:

                Alignment folder: ${aln_folder}
                    Trees folder: ${tree_folder}
                        Outgroup: ${out_group}
                 Evolution model: ${evo_model_display}
                 Distance metric: ${distance_metric_display}
                          Lambda: ${lambda}
                  Cluster method: ${cluster_method_display}
                 Infer optimal k: ${infer_optimal_k}
                         Threads: ${thread}
                            Seed: ${iqtree_seed_number}
" >>${work_dir}/parameter_input

  # Specify directory of AMAS.py
  read -p $'\nSpecify the directory of AMAS.py (Default is $HOME/Software/amas): ' amas_dir
  if [[ -z "${amas_dir}" ]]; then
    amas_dir="$HOME/Software/amas"
  fi

  # Specify directory of probdist and distance type
  if [[ "${distance_metric}" == "pd" ]] || [[ "${distance_metric}" == "pd_ex" ]]; then
    read -p $'\nSpecify the directory of probdist (Default is $HOME/Software/probdist): ' probdist_dir
    if [[ -z "${probdist_dir}" ]]; then
      probdist_dir="$HOME/Software/probdist"
    fi
    read -p $'\nSpecify distance type\nHellinger(h)\nKullback-Leibler(kl)\nJensen-Shannon(js)\nfor probabilistic distance computation(Default:js): ' probdist_type
    if [[ -z "${probdist_type}" ]]; then
      probdist_type="js"
    fi
  fi

  # Print time stamp of beginning the analysis
  time_stamp=$(date)
  echo "
Analysis started at ${time_stamp}
--------------------------------------------------------------------------------"

  ###############################
  # Part 1. Building gene trees
  ###############################
  step=1
  echo "
Part ${step}. Gene trees construction"
  SECONDS=0
  if [[ -e "${work_dir}/${tree_folder}" ]]; then
    folder_name=$(echo "${tree_folder}")
  else
    folder_name="iqMF-1000UFB"
    mkdir ${work_dir}/${folder_name} ${work_dir}/${folder_name}/report ${work_dir}/${folder_name}/trivia ${work_dir}/${folder_name}/bs-file
    bs_setting="-bb 1000 -bnni -wbtl" #Perform Ultra Fast Bootstrapping with 1000 replicates.
    cd ${work_dir}/${aln_folder}/
    parallel --no-notice -j ${thread} "iqtree -s {} ${codon_model_flag} -ninit 200 -ntop 50 -nt 1 -seed ${iqtree_seed_number} ${bs_setting} -keep-ident -quiet" ::: *.phy

    # File sorting
    mv ${work_dir}/${aln_folder}/*ufboot ${work_dir}/${folder_name}/bs-file
    mv ${work_dir}/${aln_folder}/*iqtree ${work_dir}/${folder_name}/report/
    mv ${work_dir}/${aln_folder}/*treefile ${work_dir}/${folder_name}/
    mv ${work_dir}/${aln_folder}/*log ${work_dir}/${aln_folder}/*gz ${work_dir}/${aln_folder}/*bionj ${work_dir}/${aln_folder}/*mldist ${work_dir}/${aln_folder}/*contree ${work_dir}/${aln_folder}/*splits.nex ${work_dir}/${folder_name}/trivia
    if [[ "${codon_model}" == "Yes" ]]; then mv ${work_dir}/${aln_folder}/*parstree ${work_dir}/${folder_name}/trivia; fi
  fi

  # Extract the model information
  for i in ${work_dir}/${folder_name}/report/*iqtree; do
    name=$(basename ${i} | cut -d "." -f 1-2)
    model=$(cat ${i} | grep "Best-fit model according to BIC:" | cut -d " " -f 6)
    echo -e "${name}\t${model}" >>${work_dir}/${aln_folder}/iq-modelset
  done

  # Find a single model for permutation test
  rm -rf ${work_dir}/report-temp
  mkdir ${work_dir}/report-temp
  for i in ${work_dir}/${folder_name}/report/*iqtree; do
    name=$(basename ${i} | cut -d "." -f 1)
    cat ${i} | grep -A${model_number} "List of models sorted by BIC scores:" >${work_dir}/report-temp/${name}
  done
  rm -rf ${work_dir}/report-temp2
  mkdir ${work_dir}/report-temp2
  for i in ${work_dir}/report-temp/*; do
    name=$(basename ${i})
    line_number=$(cat ${i} | wc -l)
    sed -n "4,${line_number}p" ${i} >${work_dir}/report-temp2/${name}
  done
  rm -rf ${work_dir}/report-temp
  rm -rf ${work_dir}/report-model
  mkdir ${work_dir}/report-model
  for i in ${work_dir}/report-temp2/*; do
    name=$(basename ${i})
    cat ${i} | cut -c 1-17,30-39 >${work_dir}/report-model/${name}.model
  done
  rm -rf ${work_dir}/report-temp2
  cd ${work_dir}/report-model/
  Rscript <(echo "${dBICCalculator}") 2>&1 >/dev/null
  rm -rf ${work_dir}/report-model/*model
  rm -rf ${work_dir}/report-temp
  mkdir ${work_dir}/report-temp
  for i in ${work_dir}/report-model/*; do
    if [[ ! -e ${work_dir}/report-temp/00A.txt ]]; then
      cat ${i} | cut -d " " -f 1 >${work_dir}/report-temp/00A.txt
    fi
    name=$(basename ${i})
    cat ${i} | cut -d " " -f 2 >${work_dir}/report-temp/${name}
  done
  rm -rf ${work_dir}/report-model
  cd ${work_dir}/report-temp/
  paste -d " " * >${work_dir}/model.summary
  rm -rf ${work_dir}/report-temp
  cd ${work_dir}
  Rscript <(echo "${InferSingleModel}") 2>&1 >/dev/null
  rm -f ${work_dir}/model.summary
  single_model=$(cat ${work_dir}/bestmodel.txt | sed -n "1p" | cut -d " " -f 1)

  #Root the tree and extract the rooted trees into a file.
  step_index=1
  AfterTreesConstruction ${work_dir} ${folder_name} ${out_group} ${step_index}
  rm -rf ${work_dir}/Trees
  mkdir ${work_dir}/Trees
  mv ${work_dir}/${folder_name}/ ${work_dir}/${folder_name}_rooted/ ${work_dir}/Trees/
  cd ${work_dir}
  ((step++))
  duration=$SECONDS
  time_stamp=$(date)
  echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${time_stamp}

--------------------------------------------------------------------------------"

  ############################################################
  # Part 2. Generate permutation data/concatenation sequence
  ############################################################
  if [[ "${infer_optimal_k}" == "gori" ]]; then
    echo -e "\nPart ${step}. Generate permutation data"
    SECONDS=0
    python ${amas_dir}/AMAS.py concat -i ${aln_folder}/*phy -f phylip -d dna -u phylip -t ${work_dir}/aln_concatenation.phy >/dev/null
    rm -f ${work_dir}/partitions.txt
    PhylipFormat ${work_dir}
    seed_number=2
    until [[ $((seed_number % 2)) -gt 0 ]]; do
      seed_number=${RANDOM}
    done
    if [[ "${codon_model}" == "Yes" ]]; then
      if [[ "${permutation_method}" == "f" ]]; then
        rm -rf ${work_dir}/permutation_data
        mkdir ${work_dir}/permutation_data
        cp -p ${work_dir}/aln_concatenation.phy ${work_dir}/permutation_data/
        cd ${work_dir}/permutation_data/
        codon_permute_dir=$(pwd)
        file_name="aln_concatenation.phy"
        # Parse the file
        line_number=$(cat ${codon_permute_dir}/${file_name} | wc -l)
        info=$(cat ${codon_permute_dir}/${file_name} | sed -n "1p")
        chromosome_type=$(echo "${info}" | cut -c 1)
        if [[ "${chromosome_type}" == " " ]]; then
          len=$(echo "${info}" | cut -d " " -f 3)
        else
          len=$(echo "${info}" | cut -d " " -f 2)
        fi
        cat ${codon_permute_dir}/${file_name} | sed -n "1p" >${codon_permute_dir}/info
        seq_ending=$((len + 10))
        cat ${codon_permute_dir}/${file_name} | sed -n "2,${line_number}p" | cut -c 1-10 >${codon_permute_dir}/taxa
        cat ${codon_permute_dir}/${file_name} | sed -n "2,${line_number}p" | cut -c 11-${seq_ending} >${codon_permute_dir}/aln.ori
        # Permute codon
        i=1
        while [[ "${i}" -le 100 ]]; do
          echo "${i}" >>${codon_permute_dir}/list
          ((i++))
        done
        export -f CodonJackniffPermutation
        parallel --no-notice -j 4 CodonJackniffPermutation ::: ${codon_permute_dir} ::: ${len} :::: ${codon_permute_dir}/list
        rm -rf ${codon_permute_dir}/list ${codon_permute_dir}/taxa ${codon_permute_dir}/info ${codon_permute_dir}/aln.ori
      else
        read -p "Use existant permutation datasets?(Y/n) " usept
        if [[ "${usept}" == "Y" ]] || [[ "${usept}" == "y" ]]; then
          echo -e "\nPlease put them in the folder ${work_dir}/permutation_data"
          read -p $"Continue?(Y/n) " continue_Yes_No
          if [[ "${continue_Yes_No}" != "Y" ]] && [[ "${continue_Yes_No}" != "y" ]]; then exit; fi
        else
          echo -e "\nWARNING:The datasets generation process is expensive!"
          rm -rf ${work_dir}/permutation_data
          mkdir ${work_dir}/permutation_data
          cp -p ${work_dir}/aln_concatenation.phy ${work_dir}/permutation_data/
          cd ${work_dir}/permutation_data/
          codon_permute_dir=$(pwd)
          file_name="aln_concatenation.phy"
          # Parse the file
          line_number=$(cat ${codon_permute_dir}/${file_name} | wc -l)
          info=$(cat ${codon_permute_dir}/${file_name} | sed -n "1p")
          len=$(echo "${info}" | cut -d " " -f 2)
          cat ${codon_permute_dir}/${file_name} | sed -n "1p" >${codon_permute_dir}/info
          seq_ending=$((len + 10))
          cat ${codon_permute_dir}/${file_name} | sed -n "2,${line_number}p" | cut -c 1-10 >${codon_permute_dir}/taxa
          cat ${codon_permute_dir}/${file_name} | sed -n "2,${line_number}p" | cut -c 11-${seq_ending} >${codon_permute_dir}/aln.ori
          # Divide into codon
          rm -rf ${codon_permute_dir}/codon
          mkdir ${codon_permute_dir}/codon
          iteration=$((len / 3))
          i=1
          pos_array=""
          while [[ "${i}" -le "${iteration}" ]]; do
            codon_ending_pos=$((i * 3))
            codon_starting_pos=$((codon_ending_pos - 2))
            cat ${codon_permute_dir}/aln.ori | cut -c ${codon_starting_pos}-${codon_ending_pos} >${codon_permute_dir}/codon/${i}.pos
            pos_array=$(echo "${pos_array} ${i}.pos")
            ((i++))
          done
          # Permute the codon
          echo "${pos_array}" >${codon_permute_dir}/pos.array
          ArrayPermutation_temp=$(echo "${ArrayPermutation}" | sed "s/seednumber/${seed_number}/")
          echo -e "The seed for the generation of the 100-replicate permutation data is: ${seed_number}." >${work_dir}/impMsg.txt
          Rscript <(echo "${ArrayPermutation_temp}")
          cat ${codon_permute_dir}/*shuffle >${codon_permute_dir}/pos.result
          rm -f ${codon_permute_dir}/*shuffle
          res_line_number=$(cat ${codon_permute_dir}/pos.result | wc -l)
          n=$(cat ${codon_permute_dir}/pos.result | sed -n "1p" | grep "pos" | wc -w)
          cd ${codon_permute_dir}/codon
          l=1
          while [[ "${l}" -le "${res_line_number}" ]]; do
            pos_file1=$(cat ${codon_permute_dir}/pos.result | sed -n "${l}p" | cut -d " " -f 1)
            pos_file2=$(cat ${codon_permute_dir}/pos.result | sed -n "${l}p" | cut -d " " -f 2)
            paste -d "" ${pos_file1} ${pos_file2} >${codon_permute_dir}/codon/temp0
            w=3
            while [[ "${w}" -le "${n}" ]]; do
              pos_file=$(cat ${codon_permute_dir}/pos.result | sed -n "${l}p" | cut -d " " -f ${w})
              paste -d "" ${codon_permute_dir}/codon/temp0 ${pos_file} >${codon_permute_dir}/codon/temp
              rm -f ${codon_permute_dir}/codon/temp0
              mv ${codon_permute_dir}/codon/temp ${codon_permute_dir}/codon/temp0
              ((w++))
            done
            mv ${codon_permute_dir}/codon/temp0 ${codon_permute_dir}/temp0
            paste -d "" ${codon_permute_dir}/taxa ${codon_permute_dir}/temp0 >${codon_permute_dir}/d${l}.part
            cat ${codon_permute_dir}/info ${codon_permute_dir}/d${l}.part >${codon_permute_dir}/d${l}.phy
            rm -f ${codon_permute_dir}/d${l}.part ${codon_permute_dir}/temp0
            ((l++))
          done
          rm -rf ${codon_permute_dir}/codon/ ${codon_permute_dir}/info ${codon_permute_dir}/taxa ${codon_permute_dir}/pos* ${codon_permute_dir}/aln.ori ${codon_permute_dir}/${file_name}
        fi
      fi
    else
      fseqboot -sequence ${work_dir}/aln_concatenation.phy -outfile outfile -test o -seqtype d -rewriteformat p -reps 101 -seed ${seed_number} >/dev/null
      OTU_number=6
      MakeDatasets ${work_dir} ${OTU_number} ${seed_number}
      rm -f ${work_dir}/outfile
    fi
  else
    echo -e "\nPart ${step}. Generate concatenation sequence"
    SECONDS=0
    python ${amas_dir}/AMAS.py concat -i ${aln_folder}/*phy -f phylip -d dna -u phylip -t ${work_dir}/aln_concatenation.phy >/dev/null
    rm -f ${work_dir}/partitions.txt
    PhylipFormat ${work_dir}
  fi
  ((step++))
  duration=$SECONDS
  time_stamp=$(date)
  echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${time_stamp}

--------------------------------------------------------------------------------"

  #########################
  # Part 3. Prepare files
  #########################
  if [[ "${infer_optimal_k}" == "gori" ]]; then
    echo -e "\nPart ${step}. File preparation for permutation test"
  else
    echo -e "\nPart ${step}. File preparation for clustering analysis"
  fi
  SECONDS=0

  # Generate gene trees collection
  folder_name="iqMF-brlen"
  if [[ -d "${work_dir}/${folder_name}" ]]; then rm -rf ${work_dir}/${folder_name}; fi
  mkdir ${work_dir}/${folder_name} ${work_dir}/${folder_name}/report ${work_dir}/${folder_name}/trivia
  bs_setting="" #Don't need nodes support method.
  cd ${work_dir}/${aln_folder}/
  parallel --no-notice -j ${thread} --colsep '\t' "iqtree -s {1} ${codon_model_flag} -m {2} -ninit 200 -ntop 50 -nt 1 -seed ${iqtree_seed_number} ${bs_setting} -keep-ident -quiet" :::: iq-modelset
  mv iq-modelset ${work_dir}/

  # File sorting
  mv ${work_dir}/${aln_folder}/*iqtree ${work_dir}/${folder_name}/report/
  mv ${work_dir}/${aln_folder}/*treefile ${work_dir}/${folder_name}/
  mv ${work_dir}/${aln_folder}/*log ${work_dir}/${aln_folder}/*gz ${work_dir}/${aln_folder}/*bionj ${work_dir}/${aln_folder}/*mldist ${work_dir}/${folder_name}/trivia
  if [[ "${codon_model}" == "Yes" ]]; then mv ${work_dir}/${aln_folder}/*parstree ${work_dir}/${folder_name}/trivia; fi

  # Root the tree and extract the rooted trees into a file.
  step_index=5
  AfterTreesConstruction ${work_dir} ${folder_name} ${out_group} ${step_index}
  mv ${work_dir}/${folder_name}/ ${work_dir}/${folder_name}_rooted/ ${work_dir}/Trees/

  # Get clusters information
  rm -rf ${work_dir}/Analysis
  mkdir ${work_dir}/Analysis
  mv ${work_dir}/gnTrees_collection.tre ${work_dir}/Analysis/
  if [[ -e ${work_dir}/impMsg.txt ]]; then mv ${work_dir}/impMsg.txt ${work_dir}/Analysis/; fi
  cd ${work_dir}/Analysis
  if [[ "${distance_metric}" == "pd" ]] || [[ "${distance_metric}" == "pd_ex" ]]; then
    exact_distance=$(echo "${distance_metric}" | cut -d "_" -f 2)
    if [[ "${exact_distance}" == "ex" ]]; then
      ProbabilisticDistanceCalculator ${work_dir} ${probdist_type} ${exact_distance}
    else
      # Generate model parameters file
      for i in ${work_dir}/${aln_folder}/*; do
        name=$(basename ${i})
        file_name=$(ls ${work_dir}/Trees/iqMF-1000UFB/report | grep "${name}")
        piA=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "pi(A)" | cut -d "=" -f 2 | cut -d " " -f 2)
        if [[ -z "${piA}" ]]; then
          piA="0.25"
        fi
        piC=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "pi(C)" | cut -d "=" -f 2 | cut -d " " -f 2)
        if [[ -z "${piC}" ]]; then
          piC="0.25"
        fi
        piG=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "pi(G)" | cut -d "=" -f 2 | cut -d " " -f 2)
        if [[ -z "${piG}" ]]; then
          piG="0.25"
        fi
        piT=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "pi(T)" | cut -d "=" -f 2 | cut -d " " -f 2)
        if [[ -z "${piT}" ]]; then
          piT="0.25"
        fi
        rhoAC=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "A-C" | cut -d ":" -f 2 | cut -d " " -f 2)
        rhoAG=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "A-G" | cut -d ":" -f 2 | cut -d " " -f 2)
        rhoAT=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "A-T" | cut -d ":" -f 2 | cut -d " " -f 2)
        rhoCG=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "C-G" | cut -d ":" -f 2 | cut -d " " -f 2)
        rhoCT=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "C-T" | cut -d ":" -f 2 | cut -d " " -f 2)

        alpha=$(cat ${work_dir}/Trees/iqMF-1000UFB/report/${file_name} | grep "alpha" | cut -d ":" -f 2 | cut -d " " -f 2)
        if [[ -z "${alpha}" ]]; then
          alpha=" "
        fi
        echo "${piA} ${piC} ${piG} ${piT} ${rhoAC} ${rhoAG} ${rhoAT} ${rhoCG} ${rhoCT} ${alpha}" >>${work_dir}/Analysis/modelParamFile.txt
      done
      ProbabilisticDistanceCalculator ${work_dir} ${probdist_type} ${exact_distance}
    fi
  elif [[ "${distance_metric}" == "gd" ]]; then
    Clustering_temp=$(echo "${Clustering}" | sed "s/args1/Geodesic_distance/g" | sed "s/args3/${cluster_method}/")
    Rscript <(echo "${Clustering_temp}") 2>&1 >/dev/null
  elif [[ "${distance_metric}" == "kc" ]]; then
    Clustering_temp=$(echo "${Clustering}" | sed "s/args1/KC_distance/g" | sed "s/args2/${lambda}/" | sed "s/args3/${cluster_method}/")
    Rscript <(echo "${Clustering_temp}") 2>&1 >/dev/null
  fi

  # Summarize the cluster information for each gene in each number of cluster
  cluster_info_file_number=$(ls ${work_dir}/Analysis/*cl | wc -l)
  i=1
  header=""
  while [[ "${i}" -le "${cluster_info_file_number}" ]]; do
    header=$(echo "${header} k${i}")
    ((i++))
  done
  echo "genes${header}" >${work_dir}/Analysis/cluster.summary
  while read line; do
    i=2
    cluster_sorting=""
    while [[ "${i}" -le "${cluster_info_file_number}" ]]; do
      if [[ "${i}" -lt 10 ]]; then i_mod=$(echo "k0${i}.cl"); else i_mod=$(echo "k${i}.cl"); fi
      name=$(echo "${line}" | cut -d " " -f 1)
      cluster_index=$(cat ${work_dir}/Analysis/${i_mod} | grep "${name}" | cut -d " " -f 2)
      cluster_sorting=$(echo "${cluster_sorting} ${cluster_index}")
      ((i++))
    done
    echo "${line}${cluster_sorting}" >>${work_dir}/Analysis/cluster.summary
  done <${work_dir}/Analysis/k01.cl

  # Generate order file to guide the summary
  if [[ -e ${work_dir}/Analysis/dend.order ]]; then
    have_sexchr_data_Yes_No="n"
    line_number=$(cat ${work_dir}/Analysis/dend.order | wc -l)
    i=1
    while [[ "${i}" -le "${line_number}" ]]; do
      name=$(cat ${work_dir}/Analysis/dend.order | sed -n "${i}p")
      ns_mark=$(echo "${name}" | cut -c 1-2)
      if [[ "${ns_mark}" == "ns" ]]; then
        prefix=$(echo "${name}" | cut -d "_" -f 1)
        have_sexchr_data_Yes_No="y"
        echo "${name}" >>${work_dir}/Analysis/Neo.sexpre
      else
        echo "${name}" >>${work_dir}/Analysis/A.auto
      fi
      ((i++))
    done
    rm -f dend.order
    if [[ "${have_sexchr_data_Yes_No}" == "y" ]]; then
      cat ${work_dir}/Analysis/Neo.sexpre | sort >${work_dir}/Analysis/Neo.sex
      rm -f ${work_dir}/Analysis/Neo.sexpre
      cat ${work_dir}/Analysis/A.auto ${work_dir}/Analysis/Neo.sex >${work_dir}/Analysis/order
    else
      mv ${work_dir}/Analysis/A.auto ${work_dir}/Analysis/order
    fi
    rm -f ${work_dir}/Analysis/A.auto ${work_dir}/Analysis/Neo.sex
    if [[ ! -d clinfo ]]; then mkdir clinfo; fi
    mv *cl order clinfo/
  fi
  if [[ ! -d clinfo ]]; then mkdir clinfo && mv *cl clinfo/; fi
  if [[ "${infer_optimal_k}" == "gori" ]]; then
    rm -rf permutation_test
    mkdir permutation_test
    mv clinfo/ permutation_test/

    # Get partition file
    GeneratePartitionFile ${work_dir} ${aln_folder} permutation_test

    # Manoevour the sequence files
    cp -p ${work_dir}/aln_concatenation.phy ${work_dir}/Analysis/permutation_test/
    mv ${work_dir}/permutation_data/*phy ${work_dir}/Analysis/permutation_test/
    rmdir ${work_dir}/permutation_data
  else
    rm -rf optimal_k
    mkdir optimal_k

    # Get partition file
    GeneratePartitionFile ${work_dir} ${aln_folder} optimal_k

    # Manoevour the files
    mv ${work_dir}/aln_concatenation.phy ${work_dir}/Analysis/clinfo/ ${work_dir}/Analysis/cluster.summary ${work_dir}/Analysis/optimal_k/
  fi
  ((step++))
  duration=$SECONDS
  time_stamp=$(date)
  echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${time_stamp}

--------------------------------------------------------------------------------"
  ################################################
  # Part 4. Permutation test/Clustering analysis
  ################################################
  if [[ "${infer_optimal_k}" == "gori" ]]; then
    echo -e "\nPart ${step}. Permutation test"
    SECONDS=0
    cd ${work_dir}/Analysis/permutation_test/
    current_dir=$(pwd)

    # Create partition-sorted file
    cluster_file_number=0
    rm -rf cl_sort
    mkdir cl_sort
    for file in clinfo/*cl; do
      PartitionSorting ${file}
      file_name=$(basename ${file})
      ((cluster_file_number++))
      cluster_file[$cluster_file_number]=$(echo "${file_name}_output")
      mv ./clinfo/*cl_output ./cl_sort/
    done

    # Execute the sequence splitting
    ni=0
    for file in *phy; do
      ((ni++))
      name[$ni]=$(echo "${file}" | cut -d "." -f 1)
      mkdir ${name[$ni]}

      exei=1
      while [[ "${exei}" -le "${cluster_file_number}" ]]; do
        if [[ "${codon_model}" == "Yes" ]]; then
          python ${amas_dir}/AMAS.py split -f phylip -d dna -i ${file} -l ./cl_sort/${cluster_file[$exei]} -u phylip >/dev/null
        else
          python ${amas_dir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ./cl_sort/${cluster_file[$exei]} -u phylip >/dev/null
        fi
        k_number=$(echo "${cluster_file[$exei]}" | cut -d "." -f 1)
        mkdir ./${name[$ni]}/${k_number}
        mv ./*cluster* ./${name[$ni]}/${k_number}/
        ((exei++))
      done
      rm -f ${file}
    done

    read -p $'\nContinue to Tree Construction? (y/n): ' construct_tree_Yes_No
    #construct_tree_Yes_No="y"
    if [[ "${construct_tree_Yes_No}" == "y" ]] || [[ "${construct_tree_Yes_No}" == "Y" ]]; then
      tri=1
      while [[ "${tri}" -le "${ni}" ]]; do
        trj=1
        while [[ "${trj}" -le "${cluster_file_number}" ]]; do
          k_number=$(echo "${cluster_file[$trj]}" | cut -d "." -f 1)
          cd ${current_dir}/${name[$tri]}/${k_number}/
          for seq_file in *; do
            directory=$(pwd)
            echo "${directory}/${seq_file}" >>${work_dir}/permtdata
            echo "${directory}" >>${work_dir}/permtdir-temp
          done
          ((trj++))
        done
        ((tri++))
      done
      cd ${work_dir}
      parallel --no-notice -j ${thread} "iqtree -s {} ${codon_model_flag} -m ${single_model} -ninit 200 -ntop 50 -nt 1 -seed ${iqtree_seed_number} -keep-ident -quiet" :::: permtdata
      rm -f permtdata
      cat ${work_dir}/permtdir-temp | sort -n | uniq >${work_dir}/permtdir
      rm -f ${work_dir}/permtdir-temp
      export -f loglkhExtraction
      parallel --no-notice -j ${thread} loglkhExtraction ::: ${current_dir} :::: permtdir >/dev/null
      rm -f permtdir
      if [[ ! -d "${current_dir}/logLikelihood" ]]; then mkdir ${current_dir}/logLikelihood; fi
      mv ${current_dir}/*logLhcollection ${current_dir}/logLikelihood/
      rm -rf ${current_dir}/d*

      # Extract log likelihood values from numerous files:
      mkdir ${current_dir}/value
      cd ${current_dir}/logLikelihood
      for file in *; do
        name=$(echo "${file}")
        cat ${file} | cut -d $'\t' -f 2 >${current_dir}/value/${name}
      done
      mv ${current_dir}/value/aln_concatenation_logLhcollection ${current_dir}/value/aln_concatenation_logLkhcollection

      # Generate logLkhforR.csv file
      cd ${current_dir}/value
      Rscript <(echo "${loglkhImprovementMake}") 2>&1 >/dev/null

      # File sorting
      mv ${current_dir}/value ${current_dir}/logLikelihood
      mv ${current_dir}/logLikelihood/value/logLkhcollection.csv ${current_dir}/logLikelihood/value/logLkhforR.csv ${current_dir}/
      rm -rf ${current_dir}/logLikelihood

      # Plot graph
      cd ${current_dir}
      Rscript <(echo "${loglkhImprovementPlot}") 2>&1 >/dev/null
      cd ${work_dir}/Analysis/
      Rscript -e 'library("clustree"); dd <- read.table("cluster.summary", header = TRUE); p <- clustree(dd, prefix = "k"); pdf(NULL); ggsave("cluster_summary.svg", width = 12, height = 12)' 2>&1 >/dev/null

      # File sorting
      rm -rf ${current_dir}/a.csv ${current_dir}/amd.csv ${current_dir}/b.csv ${current_dir}/bmd.csv ${current_dir}/cl_sort
    fi
  else
    echo -e "\nPart ${step}. Infer optimal number of cluster"
    SECONDS=0
    cd ${work_dir}/Analysis/optimal_k/
    current_dir=$(pwd)
    Rscript -e 'library("clustree"); dd <- read.table("cluster.summary", header = TRUE); p <- clustree(dd, prefix = "k"); pdf(NULL); ggsave("cluster_summary.svg", width = 12, height = 12)' 2>&1 >/dev/null
    if [[ "${infer_optimal_k}" == "clustree" ]]; then
      # Visualize the cluster
      decision_index=1
      while [[ "${decision_index}" -le 1 ]]; do
        if [[ -e ${work_dir}/Analysis/cluster.visualization ]]; then rm -f ${work_dir}/Analysis/cluster.visualization; fi
        read -p $'\nCheck the cluster_summary.svg and specify the number of cluser to visualize: ' k_visualization
        if [[ "${k_visualization}" -lt 10 ]]; then k_visualization_mod=$(echo "k0${k_visualization}.cl"); else k_visualization_mod=$(echo "k${k_visualization}/cl"); fi
        echo "genes M1 M2 cluster" >${work_dir}/Analysis/cluster.visualization
        while read line; do
          name=$(echo "${line}" | cut -d " " -f 1)
          k=$(cat ${current_dir}/clinfo/${k_visualization_mod} | grep "${name}" | cut -d " " -f 2)
          echo "${line} ${k}" >>${work_dir}/Analysis/cluster.visualization
        done <${work_dir}/Analysis/2-MDS.coord
        cd ${work_dir}/Analysis/
        Rscript -e 'library(ggplot2); dd <- read.table("cluster.visualization", header = TRUE); dd$cluster <- factor(dd$cluster); p <- ggplot(data = dd, aes(x = M1,y = M2, color = cluster)) + geom_point(pch = 1, size = 3);pdf(NULL); ggsave("cluster_visualization.svg", width = 12, height = 12)'
        read -p $'\nCheck cluster_visualization.svg.\nAre you happy with the number of cluster?(y/n) ' k_visualization_decision
        if [[ "${k_visualization_decision}" == "y" ]] || [[ "${k_visualization_decision}" == "Y" ]]; then
          ((decision_index++))
        fi
      done
    fi
  fi
  ((step++))
  duration=$SECONDS
  time_stamp=$(date)
  echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${time_stamp}
--------------------------------------------------------------------------------"
  ###################################
  # Part 5. Post analysis operation
  ###################################
  echo "
Part ${step}. Post analysis operation"
  SECONDS=0
  read -p "continue? " continue_Yes_No
  #continue_Yes_No="y"
  if [[ "${continue_Yes_No}" == "y" ]]; then
    if [[ "${infer_optimal_k}" == "clustree" ]]; then
      potential_k=${k_visualization}
      cd ${current_dir}
      # Create partition-sorted file
      rm -rf ${current_dir}/cl_sort
      mkdir ${current_dir}/cl_sort
      if [[ "${potential_k}" -gt 10 ]]; then potential_k_mod=$(echo "k${potential_k}"); else potential_k_mod=$(echo "k0${potential_k}"); fi
      PartitionSorting ${current_dir}/clinfo/${potential_k_mod}.cl
      mv ${current_dir}/clinfo/${potential_k_mod}.cl_output ${current_dir}/cl_sort/

      # Execute the sequence splitting
      for file in ${current_dir}/*phy; do
        name=$(basename "${file}" | cut -d "." -f 1)
        mkdir ${current_dir}/${name}
        python ${amas_dir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ${current_dir}/cl_sort/${potential_k_mod}.cl_output -u phylip >/dev/null
        mv ${current_dir}/*cluster*phy ${current_dir}/${name}/
      done

      # Supergene tree construction
      for i in ${current_dir}/${name}/*phy; do
        tree_name=$(basename "${i}" | cut -d "." -f 1 | cut -d "-" -f 1 | cut -d "_" -f 3)
        iqtree -s ${i} ${codon_model_flag} -m ${single_model} -ninit 200 -ntop 50 -nt 1 -seed ${iqtree_seed_number} -keep-ident -quiet
        nw_reroot ${i}.treefile ${out_group} >${current_dir}/${name}/${tree_name}.tre
      done
      for t in ${current_dir}/${name}/*tre; do
        spgn_tree_name=$(basename ${t} | cut -d "." -f 1)
        spgn_tree=$(cat ${t})
        echo "${spgn_tree_name} ${spgn_tree}" >>${work_dir}/gnTrees_collection.tre
        rm -f ${t}
      done
    fi
    if [[ -e ${current_dir}/logLkhforR.csv ]]; then
      # Find potential number of clusters
      cd ${current_dir}/
      Rscript <(echo "${InferPotentialk}") 2>&1 >/dev/null
      potential_k_summary=$(cat potent_knum)
      potential_k=$(sed -n "1p" potent_knum)
      read -p $"
Potential number of clusters are
${potential_k_summary}
Continue with ${potential_k} clusters? (y/n): " potential_k_reset_Yes_No
      #potential_k_reset_Yes_No="y"
      if [[ -z "${potential_k_reset_Yes_No}" ]] || [[ "${potential_k_reset_Yes_No}" == "Y" ]] || [[ "${potential_k_reset_Yes_No}" == "y" ]]; then
        echo -e "\nMsg:The analysis continues with ${potential_k} clusters"
      elif [[ "${potential_k_reset_Yes_No}" == "n" ]] || [[ "${potential_k_reset_Yes_No}" == "N" ]]; then
        read -p "Reset the number of clusters to: " potential_k_reset
        if [[ -n "${potential_k_reset}" ]] && [[ "${potential_k_reset}" -eq "${potential_k_reset}" ]] 2>/dev/null && [[ "${potential_k_reset}" -le "${cluster_file_number}" ]]; then
          potential_k=${potential_k_reset}
          echo -e "\nMsg:Number of clusters has been reset to ${potential_k} and the analysis hence continues."
        else
          echo -e "\nMsg:FAILED to reset.It continues with ${potential_k} clusters"
        fi
      else
        echo -e "\nMsg:FAILED to reset.It continues with ${potential_k} clusters"
      fi
      if [[ -n ${potential_k} ]]; then
        if [[ "${potential_k}" -gt 10 ]]; then
          potential_k_mod=$(echo "k${potential_k}")
        else
          potential_k_mod=$(echo "k0${potential_k}")
        fi
        rm -rf ${work_dir}/supergn_tree
        mkdir ${work_dir}/supergn_tree
        mv ${current_dir}/aln_concatenation/${potential_k_mod}/*treefile ${work_dir}/supergn_tree
        step_index=9
        # Extract the supergn_tree
        AfterTreesConstruction ${work_dir} supergn_tree ${out_group} ${step_index}
        mv ${work_dir}/supergn_tree* ${work_dir}/Trees/
      fi
    fi
    ((step++))
    duration=$SECONDS
    time_stamp=$(date)
    echo "
$(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
Time stamp: ${time_stamp}

--------------------------------------------------------------------------------"
  fi

  ########################
  # Part 6. File sorting
  ########################
  mv ${work_dir}/gnTrees_collection.tre ${work_dir}/supergnTrees_collection.tre
  mv ${work_dir}/supergnTrees_collection.tre ${work_dir}/Analysis/gnTrees_collection.tre ${work_dir}/Trees/
  time_stamp=$(date)
  echo "
Congratulations!!
Analysis finished at ${time_stamp}"
else
  echo -e "\nNo analysis has been run. Goodbye!"
fi
