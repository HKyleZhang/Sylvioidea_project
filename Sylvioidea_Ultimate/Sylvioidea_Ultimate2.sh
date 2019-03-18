#!/bin/bash

# Function 1. Manual
usage() {
  echo "Usage: [-a] alignment folder
       [-b] number of bootstrapping datasets
       [-c] substitution model: "nu\(nucleotide models\)", "cd\(codon models\)", "lm\(nucleotide with additional Lie Markov models\)"
       [-d] distance metirc: "gd\(Geodesic\)", "kc\(Kendall-Colijn\)", "pd\(Probabilistic distance\)"
       [-h] show usage
       [-L] a lambda between 0 and 1; NA in Geodesic and Probabilistic distance
       [-m] cluster method: "ward\(Default\)", "kmeans", "MDSK", "pam", "MDSP", "NMDSK"
       [-k] infer optimal number of clusters: "gori\(Default\)", "clustree"
       [-o] specify outgroup
       [-p] "f" to use fast method; "s" to use slow method; only valid if [-c cd]
       [-t] trees folder
       [-s] seed number for all tree construction
       [-x] number of threads"
}

# Function 2. A Fast method to permute based on codon
CodonFastPermutation() {
  codon_permute_dir=${1}
  seq_num=${2}
  seq_len=${3}
  d=${4}
  codon_sum=$((seq_len / 3))
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

  i=1
  while [[ "${i}" -le "${seq_num}" ]]; do
    label=$(cat "${codon_permute_dir}/taxa" | sed -n "${i}p")
    seq=$(cat ${codon_permute_dir}/aln.ori-${d} | sed -n "${i}p")
    echo -e "${label}\n${seq}" >>${codon_permute_dir}/d${d}.fasta
    ((i++))
  done
  rm -f ${codon_permute_dir}/aln.ori-${d}
}

# Function 3. Probabilistic distance
ProbabilisticDistanceCalculator() {
  dir=${1}
  distance_type=$(echo "-${2}")
  exact_distance=${3}

  ## Approximation type
  if [[ "${2}" == "h" ]]; then
    approx_type="-ecv"
  else
    approx_type="-e"
  fi

  ## Tree file modification
  cat ${dir}/Analysis/IQ-trees.nw | cut -d " " -f 1 >tree_name.list
  cat ${dir}/Analysis/IQ-trees.nw | cut -d " " -f 2 >gnTrees_collection_mod.tre

  ## Generate a file supplement the self-comparison distance
  tree_number=$(cat ${dir}/Analysis/IQ-trees.nw | wc -l)
  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    echo -e "\nTree ${i}, Tree ${i}, Alignment length NA, Distance 0\n" >>${dir}/Analysis/sc_dist.prob
    ((i++))
  done

  ## Calculate Probabilistic Distance
  redo_calculation=0
  while [[ "${redo_calculation}" -lt 1 ]]; do
    if [[ "${distance_type}" == "-h" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Hellinger distance." >>${work_dir}/parameter_input
    elif [[ "${distance_type}" == "-kl" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Kullback-Leibler divergence." >>${work_dir}/parameter_input
    elif [[ "${distance_type}" == "-js" ]]; then
      echo -e "\nNB:Probabilistic distance is calcalated in Jensen-Shannon distance." >>${work_dir}/parameter_input
    fi

    ## Calculation begins
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

  ## File modification
  cat tree_dist_raw.prob sc_dist.prob >tree_dist.prob
  rm -f tree_dist_raw.prob sc_dist.prob gnTrees_collection_mod.tre
  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    i_mod=$(printf "%0${#tree_number}d" "${i}")
    sed -i "s/Tree ${i},/Tree${i_mod}/g" tree_dist.prob
    ((i++))
  done
  Rscript <(echo "${ProbDistMatrix}") 2>&1 >/dev/null

  i=1
  while [[ "${i}" -le "${tree_number}" ]]; do
    i_mod=$(printf "%0${#tree_number}d" "${i}")
    tree_name=$(cat ${dir}/Analysis/tree_name.list | sed -n "${i}p")
    sed -i "s/${i_mod}/${tree_name}/g" ${dir}/Analysis/dist_matrix.txt
    ((i++))
  done
  rm -f tree_name.list

  ## Clustering
  Rscript <(echo "${ClusteringProbdist}") 2>&1 >/dev/null
}

# Function 4. Generate file containing partitions
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

# Function 5. Sorting out the partitions according to the clustering results
PartitionSorting() {
  output_name=$(echo "${2}_output")
  line_number=$(cat ${1}/partition | wc -l)

  i=1
  while [[ "${i}" -le "${line_number}" ]]; do
    gene_name[$i]=$(cat ${1}/partition | sed -n "${i}p" | cut -d " " -f 1)
    partition[$i]=$(cat ${1}/partition | sed -n "${i}p" | cut -d " " -f 3 | cut -d ";" -f 1)
    ((i++))
  done

  i=1
  while [[ "${i}" -le "${line_number}" ]]; do
    cluster[$i]=$(cat ${1}/clinfo/${2} | grep "${gene_name[$i]}" | cut -d " " -f 2)
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

      k_mod=$(printf "%02d" "${k}")
      echo "DNA, cluster${k_mod} = ${partition_content}" >>${1}/clinfo/${output_name}
      ((k++))
    fi

    ((i++))
  done
}

# R code 1. Infer the single model
InferSingleModel='
suppressPackageStartupMessages(library(tidyverse))

diff_sq <- function(df, min_bic)
  mutate_at(df, vars(dBIC), list(~(. - min_bic) ^ 2))

model_files <- list.files(pattern = "*.model")
col_names <- c("Model", "dBIC")
model_sheets <-
  map(model_files, read.table) %>%
  map(setNames, col_names) %>%
  set_names(nm = model_files)

BIC_min <- map(model_sheets, function(df)
  min(df$dBIC))
model_sheets <-
  suppressWarnings(map2(model_sheets, BIC_min, diff_sq) %>%
                     reduce(full_join, by = "Model")) %>%
  column_to_rownames(var = "Model") %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Model") %>%
  setNames(col_names)
model_sheets[which.min(model_sheets$dBIC), 1]'

# R code 2. Alignment concatenation
aln_concat='
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ape))

aln_file <- list.files(pattern = "*.phy")
phylip_aln <- map(aln_file, read.dna, format = "sequential")
concat <- phylip_aln %>% reduce(cbind.DNAbin)
write.FASTA(concat, "aln_concatenation.fasta")'

# R code 3. Permute array
ArrayPermutation='
array <- read.table("pos.array")
seed_number <- SEEDNUM
set.seed(seed_number)
for (i in 1:BSNUM) {
  seed_number <- seed_number + 1
  set.seed(seed_number)
  p <- sample(array)
  name <- paste(i, ".shuffle", sep = "")
  write.table(
    p,
    file = name,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}'

# R code 4. Root phylogenetic trees with phytools R package
ReRooTree='
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyverse))

reroot_mod <- function(t) {
  outgroup_index <- grep("OUTGROUP", t$tip.label)
  reroot(t, outgroup_index)
}

iqtree_file <- list.files(pattern = "*.treefile")
iqtree_label <-
  list(iqtree_file) %>% map(gsub, pattern = ".phy.treefile", replacement = " ") %>% unlist()
iqtrees <-
  lapply(iqtree_file, read.tree) %>% set_names(iqtree_label) %>% lapply(reroot_mod)
class(iqtrees) <- "multiPhylo"
write.tree(iqtrees, "IQ-trees.nw", tree.names = TRUE)'

# R code 5. Root phylogenetic trees with ape R package
RooTree='
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(tidyverse))

iqtree_file <- list.files(pattern = "*.treefile")
iqtree_label <-
  list(iqtree_file) %>% map(gsub, pattern = ".phy.treefile", replacement = "") %>% unlist()
iqtrees <-
  lapply(iqtree_file, read.tree) %>% set_names(iqtree_label)
class(iqtrees) <- "multiPhylo"
iqtrees <- root.multiPhylo(iqtrees, "OUTGROUP", resolve.root = T)
write.tree(iqtrees, "IQ-trees.nw", tree.names = TRUE)'

# R code 6. Make ProbDist matrix
ProbDistMatrix='
suppressPackageStartupMessages(library(tidyverse))

dd <- read.table("tree_dist.prob")
dd <- data.frame(dd$V1, dd$V2, dd$V7)
names(dd) <- c("a", "b", "value")
dm <-
  spread(dd, a, value) %>% 
  column_to_rownames(var = "b") %>% 
  as.dist() %>% 
  as.matrix()
write.table(dm, "dist_matrix.txt", quote = F)'

# R code 7. Get clusters
Clustering='
suppressMessages(library(treespace))

trees <- read.tree("IQ-trees.nw")
cluster_number <- length(trees)

if ("args1" == "KC_distance") {
  L <- args2
  KC_res <-
    treespace(
      trees,
      method = "treeVec",
      nf = 2,
      lambda = L,
      return.tree.vectors = TRUE
    )
  distance <- KC_res$D
} else if ("args1" == "Geodesic_distance") {
  suppressPackageStartupMessages(library(distory))
  suppressPackageStartupMessages(library(phangorn))

  geod <- dist.multiPhylo(trees)
  wrfd <- wRF.dist(trees)
  
  wrfd <- as.matrix(wrfd)
  treename <- rownames(wrfd)
  geod <- as.matrix(geod)
  rownames(geod) <- treename
  colnames(geod) <- treename
  distance <- as.dist(geod)
}

if ("args3" == "ward") {
  ## Hierachichal clustering with MDS visualization
  distance.hclust <- hclust(distance, method = "ward.D2")
  svg("dendrogram.svg", width = 10, height = 10)
  plot(distance.hclust, cex = 0.6)
  dev.off()
  dend_label_raw <- as.data.frame(distance.hclust$labels)
  dend_label <-
    dend_label_raw$`distance.hclust$labels`[order.dendrogram(as.dendrogram(distance.hclust))]
  dend_label <- as.data.frame(dend_label)
  write.table(
    dend_label,
    file = "dend.order",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )

  if (cluster_number > 15) {
    cluster_number <- 15
  }

  for (j in 1:cluster_number) {
    cl <- cutree(distance.hclust, k = j)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(cl, nj, quote = FALSE, col.names = FALSE)
  }
}

## K-means on tree vectors
if (("args1" == "KC_distance") & ("args3" == "kmeans")) {
  treVec <- KC_res$vectors

  if (cluster_number > 15) {
    cluster_number <- 15
  } else {
    cluster_number <- cluster_number - 1
  }

  for (j in 1:cluster_number) {
    kc.km <- kmeans(treVec, centers = j, nstart = 50)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
  }
}

## Multi-Dimensional scaling down to 2D and 3D, K-means on 3D MDS
if ("args3" == "MDSK") {
  suppressPackageStartupMessages(library(smacof))

  mdsres2D <- smacofSym(distance, ndim = 2, type = "ratio")
  treMDS2D <- as.data.frame(mdsres2D$conf)
  write.table(treMDS2D,
              file = "2-MDS.coord",
              col.names = FALSE,
              quote = FALSE)
  write.table(
    mdsres2D$stress,
    file = "2-MDS.stress",
    col.names = FALSE,
    quote = FALSE
  )
  mdsres3D <- smacofSym(distance, ndim = 3, type = "ratio")
  treMDS3D <- as.data.frame(mdsres3D$conf)
  write.table(treMDS3D,
              file = "3-MDS.coord",
              col.names = FALSE,
              quote = FALSE)
  write.table(
    mdsres3D$stress,
    file = "3-MDS.stress",
    col.names = FALSE,
    quote = FALSE
  )

  if (cluster_number > 15) {
    cluster_number <- 15
  } else {
    cluster_number <- cluster_number - 1
  }

  for (j in 1:cluster_number) {
    kc.km <- kmeans(treMDS3D, centers = j, nstart = 50)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
  }
}

## PAM on the distance/dissimilarities matrix
if ("args3" == "pam") {
  suppressPackageStartupMessages(library(cluster))

  if (cluster_number > 15) {
    cluster_number <- 15
  } else {
    cluster_number <- cluster_number - 1
  }

  for (j in 1:cluster_number) {
    kc.pam <- pam(distance,
                  diss = TRUE,
                  k = j,
                  cluster.only = TRUE)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(kc.pam, nj, quote = FALSE, col.names = FALSE)
  }
}

# Multi-Dimensional scaling down to 2D and 3D, PAM on 3D MDS
if ("args3" == "MDSP") {
  suppressPackageStartupMessages(library(smacof))
  suppressPackageStartupMessages(library(cluster))

  mdsres2D <- smacofSym(distance, ndim = 2, type = "ratio")
  treMDS2D <- as.data.frame(mdsres2D$conf)
  write.table(treMDS2D,
              file = "2-MDS.coord",
              col.names = FALSE,
              quote = FALSE)
  write.table(
    mdsres2D$stress,
    file = "2-MDS.stress",
    col.names = FALSE,
    quote = FALSE
  )
  mdsres3D <- smacofSym(distance, ndim = 3, type = "ratio")
  treMDS3D <- as.data.frame(mdsres3D$conf)
  write.table(treMDS3D,
              file = "3-MDS.coord",
              col.names = FALSE,
              quote = FALSE)
  write.table(
    mdsres3D$stress,
    file = "3-MDS.stress",
    col.names = FALSE,
    quote = FALSE
  )

  if (cluster_number > 15) {
    cluster_number <- 15
  } else {
    cluster_number <- cluster_number - 1
  }
  for (j in 1:cluster_number) {
    kc.pam <- pam(treMDS3D,
                  diss = FALSE,
                  k = j,
                  cluster.only = TRUE)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(kc.pam, nj, quote = FALSE, col.names = FALSE)
  }
}

## Non-metric Multi-Dimensional scaling down to 2D and K-means
if ("args3" == "NMDSK") {
  suppressPackageStartupMessages(library(smacof))

  mdsres2D <- smacofSym(distance, ndim = 2, type = "ordinal")
  treMDS2D <- as.data.frame(mdsres2D$conf)
  write.table(treMDS2D,
              file = "2-MDS.coord",
              col.names = FALSE,
              quote = FALSE)
  write.table(
    mdsres2D$stress,
    file = "2-MDS.stress",
    col.names = FALSE,
    quote = FALSE
  )

  if (cluster_number > 15) {
    cluster_number <- 15
  } else {
    cluster_number <- cluster_number - 1
  }
  for (j in 1:cluster_number) {
    kc.km <- kmeans(treMDS2D, centers = j, nstart = 50)
    if (j < 10) {
      j = paste("0", j, sep = "")
      nj = paste("k", j, ".cl", sep = "")
    } else {
      nj = paste("k", j, ".cl", sep = "")
    }
    write.table(kc.km$cluster, nj, quote = FALSE, col.names = FALSE)
  }
}'

# R code 8. Graph plotting for permutation test
dloglkhPlot='
suppressPackageStartupMessages(library(tidyverse))

dloglkh <- function(df)
  mutate_at(df, vars(V3), list(~(lead(.) - .))) %>% filter(is.na(V3) == FALSE)

# Import and split loglkh.collection
dd <- read.table("loglkh.collection")
dd_actual <- dd[dd$V1 == "actual_data", ]
dd_permute <- dd[dd$V1 != "actual_data", ]

# Calculate differece in log likelihood between different number of clusters
dd_actual <-
  mutate_at(dd_actual, vars(V3), list( ~ (lead(.) - .))) %>% select(V2, V3) %>% filter(is.na(V3) == FALSE)
colnames(dd_actual) <- c("dk", "dloglkh")

dd_p_nested <- dd_permute %>% group_by(V1) %>% nest()
dd_p_nested <- dd_p_nested$data %>% map(dloglkh)

dd_p_num <- length(dd_p_nested)

dd_p_nested <-
  dd_p_nested %>% reduce(full_join, by = "V2") %>% column_to_rownames(var = "V2")
colnames(dd_p_nested) <- seq(1:dd_p_num)

dd_permute <-
  dd_p_nested %>% t() %>% as.data.frame() %>% gather(dk, dloglkh)
dd_permute_mean <-
  dd_p_nested %>% add_column(mean_dloglkh = apply(dd_p_nested, 1, mean)) %>% select(mean_dloglkh) %>% rownames_to_column(var = "dk")

# Plotting
xlabel <- "1→2"
for (i in 2:nrow(dd_actual)) {
  xlabel_temp <- paste(i, "→", i + 1, sep = "")
  xlabel <- c(xlabel, xlabel_temp)
}

p <- ggplot() +
  geom_boxplot(
    data = dd_permute,
    aes(x = dk, y = dloglkh),
    color = "grey"
  ) +
  geom_dotplot(
    data = dd_permute,
    aes(x = dk, y = dloglkh),
    binwidth = 0.1,
    binaxis = "y",
    stackdir = "center",
    fill = "grey",
  ) +
  geom_line(data = dd_actual, aes(
    x = dk,
    y = dloglkh,
    group = 1,
    color = "dk"
  )) +
  geom_point(data = dd_actual, aes(x = dk, y = dloglkh), color = "firebrick2") +
  geom_line(data = dd_permute_mean, aes(
    x = dk,
    y = mean_dloglkh,
    group = 1,
    color = "dloglkh"
  )) +
  geom_point(data = dd_permute_mean, aes(x = dk, y = mean_dloglkh), color = "royalblue2") +
  scale_x_discrete(labels = xlabel) +
  labs(x = "Number of clusters, k", y = "log Likelihood improvement") +
  theme(
    axis.text.x = element_text(
      size = 18,
      vjust = 0.5,
      angle = 45,
      family = "Calibri"
    ),
    axis.title.x = element_text(size = 20, family = "Calibri"),
    axis.text.y = element_text(size = 18, family = "Calibri"),
    axis.title.y = element_text(size = 20, family = "Calibri"),
    legend.title = element_text(size = 20, family = "Calibri"),
    legend.text = element_text(size = 20, family = "Calibri")
  ) + scale_color_manual(
    name = "Legend",
    values = c("firebrick2", "royalblue2"),
    labels = c("Actual data", "Permutation data mean")
  )

pdf(NULL)
ggsave("permutation_test.svg",
       p,
       width = 14,
       height = 10
)'

# Pre-step: Load the flags and check the settings #
###################################################
while getopts "ha:b:c:d:u:L:m:k:o:p:s:t:x:" opt; do
  case $opt in
  h) usage && exit ;;
  a) aln_folder=$(echo "${OPTARG}" | cut -d "/" -f 1) ;;
  b) bs_num=${OPTARG} ;;
  c) sub_model=${OPTARG} ;;
  d) distance_metric=${OPTARG} ;;
  L) lambda=${OPTARG} ;;
  m) cluster_method=${OPTARG} ;;
  k) infer_optimal_k=${OPTARG} ;;
  o) out_group=${OPTARG} ;;
  p) permutation_method=${OPTARG} ;;
  s) iqtree_seed_number=${OPTARG} ;;
  t) tree_folder=$(echo "${OPTARG}" | cut -d "/" -f 1) ;;
  x) thread=${OPTARG} ;;
  esac
done

## No argument returns usage manual and exit ##
if [[ $OPTIND -eq 1 ]]; then
  usage
  exit
fi

## Get current directory ##
work_dir=$(pwd)

## Check if alignment folder ##
if [[ -z "${aln_folder}" ]] || [[ ! -d "${work_dir}/${aln_folder}" ]]; then
  echo -e "\nMsg: Some files are missing.\n"
  usage
  exit
fi

## Check trees folder ##
if [[ -z "${tree_folder}" ]]; then
  tree_folder="NA"
fi

## Check outgroup ##
if [[ -z "${out_group}" ]]; then
  file=$(ls ${work_dir}/${aln_folder} | sed -n "1p")
  out_group=$(cat ${work_dir}/${aln_folder}/${file} | sed -n "2p" | cut -d " " -f 1)
fi

## Check distance metric ##
if [[ -z "${distance_metric}" ]]; then
  distance_metric="kc"
fi

## Check substitution model ##
if [[ "${sub_model}" == "nu" ]] || [[ -z "${sub_model}" ]]; then
  model_number=134
  add_model_flag=""
  evo_model_display="Nucleotide substitution models"
elif [[ "${sub_model}" == "cd" ]]; then
  model_number=92
  add_model_flag="-st CODON"
  evo_model_display="Codon substitution models"
  if [[ "${distance_metric}" == "pd" ]]; then
    echo -e "\nMsg: Probabilistic distance is NOT applicable when codon model is chosen."
    exit
  fi
elif [[ "${sub_model}" == "lm" ]]; then
  model_number=358
  add_model_flag="-m MFP+LMRY"
  evo_model_display="Lie Markov models"
fi

## Check if probabilistic distance is used ##
if [[ "${distance_metric}" == "pd" ]] || [[ "${distance_metric}" == "pd_ex" ]]; then
  distance_metric_display="Probabilistic distance"
  lambda="NA"
elif [[ "${distance_metric}" == "gd" ]]; then
  distance_metric_display="Geodesic distance"
  lambda="NA"
elif [[ "${distance_metric}" == "kc" ]]; then
  distance_metric_display="Kendall-Colijn metric"
fi

## Check the lambda setting ##
if [[ "${distance_metric}" == "kc" ]]; then
  if [[ -z "${lambda}" ]]; then
    echo -e "\nERROR: Missing lambda value!"
    exit
  fi
fi

## Check the clustering method setting ##
if [[ -z "${cluster_method}" ]]; then
  cluster_method="ward"
fi

## Cluster method info ##
case ${cluster_method} in
ward) cluster_method_display="Hierachichal clustering" ;;
kmeans) cluster_method_display="k-means clustering" ;;
MDSK) cluster_method_display="MDS -> k-means clustering" ;;
pam) cluster_method_display="k-medoids (PAM) clustering" ;;
MDSP) cluster_method_display="MDS -> k-medoids (PAM) clustering" ;;
NMDSK) cluster_method_display="NMDS" ;;
esac

## Check method for inferring optimal cluster number ##
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

## Check number of bootstrapping datasets ##
if [[ "${infer_optimal_k}" == "clustree" ]]; then
  bs_num="NA"
else
  if [[ -z "${bs_num}" ]]; then
    bs_num=100
  fi
fi

## Check number of threads ##
if [[ -z "${thread}" ]]; then
  thread=$(nproc --all)
fi

## Check the seed for IQ-Tree ##
if [[ -z "${iqtree_seed_number}" ]]; then
  iqtree_seed_number=${RANDOM}
fi

## Show the requirements and settings ##
Programme_header='
=========================================================================
                        SYLVIOIDEA Ultimate
======================================== version 2 ======================'

input_setting=$(echo -e "${Programme_header}

SETTINGS:
                Alignment folder: ${aln_folder}
                    Trees folder: ${tree_folder}
                        Outgroup: ${out_group}
                 Evolution model: ${evo_model_display}
                 Distance metric: ${distance_metric_display}
                          Lambda: ${lambda}
                  Cluster method: ${cluster_method_display}
                 Infer optimal k: ${infer_optimal_k}
      No. bootstrapping datasets: ${bs_num}
                         Threads: ${thread}
                            Seed: ${iqtree_seed_number}\n")

read -p "${input_setting}
Continue?(Y/n): " answer
if [[ "${answer}" != "Y" ]] && [[ "${answer}" != "y" ]]; then
  echo -e "\nNo analysis has been run. Goodbye!"
  exit
fi

# Programme starts #
####################
echo "${input_setting}" >>${work_dir}/Input_setting

## Specify directory of probdist and distance type ##
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

time_stamp=$(date)
echo -e "\nAnalysis started at\n${time_stamp}.\n--------------------\n"

# Part 1. Phylogenetic trees construction #
###########################################
step=1
echo "Part ${step}. Phylogenetic trees construction"

if [[ -e "${work_dir}/${tree_folder}" ]]; then
  folder_name=$(echo "${tree_folder}")
else
  folder_name="IQ-Tree_output"
  mkdir ${work_dir}/${folder_name} \
    ${work_dir}/${folder_name}/report \
    ${work_dir}/${folder_name}/trivia
  cd ${work_dir}/${aln_folder}/

  if [[ "${distance_metric}" == "kc" ]]; then
    parallel --no-notice -j ${thread} "iqtree -s {} ${add_model_flag} --polytomy --merit BIC -nt 1 -seed ${iqtree_seed_number} -keep-ident -quiet" ::: *.phy
  else
    parallel --no-notice -j ${thread} "iqtree -s {} ${add_model_flag} --merit BIC -nt 1 -seed ${iqtree_seed_number} -keep-ident -quiet" ::: *.phy
  fi

  mv ${work_dir}/${aln_folder}/*iqtree ${work_dir}/${folder_name}/report/
  mv ${work_dir}/${aln_folder}/*treefile ${work_dir}/${folder_name}/
  mv ${work_dir}/${aln_folder}/*log ${work_dir}/${aln_folder}/*gz ${work_dir}/${aln_folder}/*bionj ${work_dir}/${aln_folder}/*mldist ${work_dir}/${folder_name}/trivia

  if [[ "${sub_model}" == "cd" ]]; then
    mv ${work_dir}/${aln_folder}/*parstree ${work_dir}/${folder_name}/trivia
  fi
fi

## Find a single model for permutation test ##
mkdir ${work_dir}/report-temp

for i in ${work_dir}/${folder_name}/report/*iqtree; do
  name=$(basename ${i} | cut -d "." -f 1)
  cat ${i} | grep -A${model_number} "List of models sorted by BIC scores:" | sed -n "4,${model_number}p" | awk '{print $1, $9}' >${work_dir}/report-temp/${name}.model
done

cd ${work_dir}/report-temp/
single_model=$(Rscript <(echo "${InferSingleModel}") | cut -d "\"" -f 2)
echo -e "\n  Message: \"${single_model}\" will be used for \n  tree reconstruction in permutation test."
rm -rf ${work_dir}/report-temp/

((step++))
echo -e "\n--------------------"

# Part 2. Concatenation alignment/Generate Permutation dataset #
################################################################
if [[ "${infer_optimal_k}" == "clustree" ]]; then
  echo -e "\nPart ${step}. Generate concatenation alignment"

  cd ${work_dir}/${aln_folder}
  Rscript <(echo "${aln_concat}")
  mv ${work_dir}/${aln_folder}/aln_concatenation.fasta ${work_dir}
else
  echo -e "\nPart ${step}. Generate permutation dataset"

  cd ${work_dir}/${aln_folder}
  Rscript <(echo "${aln_concat}")
  mv ${work_dir}/${aln_folder}/aln_concatenation.fasta ${work_dir}

  seed_number=2
  until [[ $((seed_number % 2)) -gt 0 ]]; do
    seed_number=${RANDOM}
  done

  if [[ "${sub_model}" == "cd" ]]; then
    mkdir ${work_dir}/Permutation_dataset
    concatenation_aln="aln_concatenation.fasta"
    cp -p ${work_dir}/${concatenation_aln} ${work_dir}/Permutation_dataset/
    codon_permute_dir=$(echo "${work_dir}/Permutation_dataset")

    ## Parse the file ##
    cat ${codon_permute_dir}/${concatenation_aln} | grep "^>" | sed "/^$/d" >${codon_permute_dir}/taxa
    cat ${codon_permute_dir}/${concatenation_aln} | sed "s/^>.*//g" | sed "/^$/d" >${codon_permute_dir}/aln.ori
    seq_num=$(cat ${codon_permute_dir}/taxa | grep -c "^>")
    seq_len=$(cat ${codon_permute_dir}/${concatenation_aln} | sed -n "2p" | wc -c)

    ## Generation of permutation datasets ##
    if [[ "${permutation_method}" == "f" ]]; then
      i=1
      while [[ "${i}" -le "${bs_num}" ]]; do
        echo "${i}" >>${codon_permute_dir}/list
        ((i++))
      done

      export -f CodonFastPermutation
      parallel --no-notice -j ${thread} CodonFastPermutation ::: ${codon_permute_dir} ::: ${seq_num} ::: ${seq_len} :::: ${codon_permute_dir}/list
      rm -rf ${codon_permute_dir}/list ${codon_permute_dir}/taxa ${codon_permute_dir}/aln.ori
    else
      read -p "Use existant permutation datasets?(Y/n) " usept
      if [[ "${usept}" == "Y" ]] || [[ "${usept}" == "y" ]]; then
        echo -e "\nPlease put them in the folder ${work_dir}/Permutation_dataset"
        read -p $"Continue?(Y/n) " continue_Yes_No
        if [[ "${continue_Yes_No}" != "Y" ]] && [[ "${continue_Yes_No}" != "y" ]]; then exit; fi
      else
        ### Divide into codons ###
        mkdir ${codon_permute_dir}/codon
        iteration=$((seq_len / 3))

        i=1
        pos_array=""
        while [[ "${i}" -le "${iteration}" ]]; do
          codon_ending_pos=$((i * 3))
          codon_starting_pos=$((codon_ending_pos - 2))
          cat ${codon_permute_dir}/aln.ori | cut -c ${codon_starting_pos}-${codon_ending_pos} >${codon_permute_dir}/codon/${i}.pos
          pos_array=$(echo "${pos_array} ${i}.pos")
          ((i++))
        done

        ### Permute the codon ###
        echo "${pos_array}" >${codon_permute_dir}/pos.array
        echo -e "The seed for the generation of the ${bs_num}-replicate permutation data is: ${seed_number}." >${work_dir}/Important_Message.txt
        cd ${codon_permute_dir}
        ArrayPermutation_temp=$(echo "${ArrayPermutation}" | sed "s/SEEDNUM/${seed_number}/" | sed "s/BSNUM/${bs_num}/")
        Rscript <(echo "${ArrayPermutation_temp}")

        cat ${codon_permute_dir}/*shuffle >${codon_permute_dir}/pos.result
        rm -f ${codon_permute_dir}/*shuffle
        res_line_number=$(cat ${codon_permute_dir}/pos.result | wc -l)
        n=$(cat ${codon_permute_dir}/pos.result | sed -n "1p" | grep -o "pos" | wc -l)
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

          i=1
          while [[ "${i}" -le "${seq_num}" ]]; do
            label=$(cat "${codon_permute_dir}/taxa" | sed -n "${i}p")
            seq=$(cat ${codon_permute_dir}/temp0 | sed -n "${i}p")
            echo -e "${label}\n${seq}" >>${codon_permute_dir}/d${l}.fasta
            ((i++))
          done

          rm -f ${codon_permute_dir}/temp0
          ((l++))
        done
        rm -rf ${codon_permute_dir}/codon/ ${codon_permute_dir}/taxa ${codon_permute_dir}/pos* ${codon_permute_dir}/aln.ori ${codon_permute_dir}/${concatenation_aln}
      fi
    fi
  else
    ### Nucleotide permutation ###
    fseqboot -sequence ${work_dir}/aln_concatenation.fasta -outfile ${work_dir}/outfile -test o -seqtype d -rewriteformat p -reps $((bs_num + 1)) -seed ${seed_number} 2>&1 >/dev/null
    echo -e "The seed for the generation of the ${bs_num}-replicate permutation data is: ${seed_number}." >${work_dir}/Important_Message.txt

    seq_num=$(cat ${work_dir}/aln_concatenation.fasta | grep -c "^>")
    file_number=1
    mkdir ${work_dir}/Permutation_dataset

    i=1
    j=$((i + bs_num))
    while [[ "${i}" -lt "${j}" ]]; do
      k=$((i + 1))
      start_number=$(cat ${work_dir}/outfile | grep -n "${seq_num}" | cut -d ":" -f 1 | sed -n "${i}p")
      end=$(cat ${work_dir}/outfile | grep -n "${seq_num}" | cut -d ":" -f 1 | sed -n "${k}p")
      end_number=$((end - 1))
      file_num_mod=$(printf "%0${#bs_num}d" "${file_number}")
      cat ${work_dir}/outfile | sed -n "${start_number},${end_number}p" >${work_dir}/Permutation_dataset/d${file_num_mod}.phy
      ((i++))
      ((file_number++))
    done

    rm -f ${work_dir}/outfile
  fi
fi
((step++))
echo -e "\n--------------------"

# Part 3. Permutation test/Clustering analysis #
################################################

if [[ "${infer_optimal_k}" == "gori" ]]; then
  echo -e "\nPart ${step}. Permutation test"
else
  echo -e "\nPart ${step}. Clustering analysis"
fi

## Get clusters information ##
mkdir ${work_dir}/Analysis

if [[ -e ${work_dir}/Important_Message.txt ]]; then
  mv ${work_dir}/Important_Message.txt ${work_dir}/Analysis/
fi

cd ${work_dir}/IQ-Tree_output
ReRooTree_temp=$(echo "${ReRooTree}" | sed "s/OUTGROUP/${out_group}/")
Rscript <(echo "${ReRooTree_temp}")
mv ${work_dir}/IQ-Tree_output/IQ-trees.nw ${work_dir}/Analysis/
cd ${work_dir}/Analysis

if [[ "${distance_metric}" == "pd" ]] || [[ "${distance_metric}" == "pd_ex" ]]; then
  exact_distance=$(echo "${distance_metric}" | cut -d "_" -f 2)

  if [[ "${exact_distance}" == "ex" ]]; then
    ProbabilisticDistanceCalculator ${work_dir} ${probdist_type} ${exact_distance}
  else
    ### Generate model parameters file for Probablistic distance ###
    for i in ${work_dir}/${aln_folder}/*; do
      name=$(basename ${i})
      file_name=$(ls ${work_dir}/IQ-Tree_output/report | grep "${name}")
      piA=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "pi(A)" | cut -d "=" -f 2 | cut -d " " -f 2)
      if [[ -z "${piA}" ]]; then
        piA="0.25"
      fi
      piC=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "pi(C)" | cut -d "=" -f 2 | cut -d " " -f 2)
      if [[ -z "${piC}" ]]; then
        piC="0.25"
      fi
      piG=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "pi(G)" | cut -d "=" -f 2 | cut -d " " -f 2)
      if [[ -z "${piG}" ]]; then
        piG="0.25"
      fi
      piT=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "pi(T)" | cut -d "=" -f 2 | cut -d " " -f 2)
      if [[ -z "${piT}" ]]; then
        piT="0.25"
      fi
      rhoAC=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "A-C" | cut -d ":" -f 2 | cut -d " " -f 2)
      rhoAG=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "A-G" | cut -d ":" -f 2 | cut -d " " -f 2)
      rhoAT=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "A-T" | cut -d ":" -f 2 | cut -d " " -f 2)
      rhoCG=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "C-G" | cut -d ":" -f 2 | cut -d " " -f 2)
      rhoCT=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "C-T" | cut -d ":" -f 2 | cut -d " " -f 2)

      alpha=$(cat ${work_dir}/IQ-Tree_output/report/${file_name} | grep "alpha" | cut -d ":" -f 2 | cut -d " " -f 2)
      if [[ -z "${alpha}" ]]; then
        alpha=" "
      fi
      echo "${piA} ${piC} ${piG} ${piT} ${rhoAC} ${rhoAG} ${rhoAT} ${rhoCG} ${rhoCT} ${alpha}" >>${work_dir}/Analysis/modelParamFile.txt
    done

    ProbabilisticDistanceCalculator ${work_dir} ${probdist_type} ${exact_distance}
  fi

elif [[ "${distance_metric}" == "gd" ]]; then
  Clustering_temp=$(echo "${Clustering}" | sed "s/args1/Geodesic_distance/g" | sed "s/args3/${cluster_method}/")
elif [[ "${distance_metric}" == "kc" ]]; then
  Clustering_temp=$(echo "${Clustering}" | sed "s/args1/KC_distance/g" | sed "s/args2/${lambda}/" | sed "s/args3/${cluster_method}/")
fi
Rscript <(echo "${Clustering_temp}")

## Summarize the cluster information for each gene in each number of cluster ##
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
    i_mod=$(printf "%02d" "${i}")
    name=$(echo "${line}" | cut -d " " -f 1)
    cluster_index=$(cat ${work_dir}/Analysis/k${i_mod}.cl | grep "${name}" | cut -d " " -f 2)
    cluster_sorting=$(echo "${cluster_sorting} ${cluster_index}")
    ((i++))
  done
  echo "${line}${cluster_sorting}" >>${work_dir}/Analysis/cluster.summary
done <${work_dir}/Analysis/k01.cl

## Generate order file to guide the summary ##
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

  if [[ ! -d ${work_dir}/Analysis/clinfo ]]; then
    mkdir ${work_dir}/Analysis/clinfo
  fi
  mv ${work_dir}/Analysis/*cl ${work_dir}/Analysis/order ${work_dir}/Analysis/clinfo/
fi

if [[ ! -d ${work_dir}/Analysis/clinfo ]]; then
  mkdir ${work_dir}/Analysis/clinfo
  mv ${work_dir}/Analysis/*cl clinfo/
fi

## Get partition file ##
if [[ "${infer_optimal_k}" == "gori" ]]; then
  mkdir ${work_dir}/Analysis/Permutation_test
  mv ${work_dir}/Analysis/clinfo/ Permutation_test/

  GeneratePartitionFile ${work_dir} ${aln_folder} Permutation_test
else
  mkdir ${work_dir}/Analysis/optimal_k

  GeneratePartitionFile ${work_dir} ${aln_folder} optimal_k

  mv ${work_dir}/Analysis/clinfo/ ${work_dir}/Analysis/cluster.summary ${work_dir}/Analysis/optimal_k/
fi

if [[ "${infer_optimal_k}" == "gori" ]]; then
  ## Permutation test starts ##
  current_dir=$(echo "${work_dir}/Analysis/Permutation_test")
  mv ${work_dir}/aln_concatenation.fasta ${work_dir}/Permutation_dataset/actual_data.fasta
  mv ${work_dir}/Permutation_dataset/ ${current_dir}

  ### Create partition-sorted file ###
  mkdir ${current_dir}/cl_sort/

  cluster_file_number=0
  for file in ${current_dir}/clinfo/*cl; do
    file_name=$(basename ${file})
    PartitionSorting ${current_dir} ${file_name}
    ((cluster_file_number++))
    cluster_file[$cluster_file_number]=$(echo "${file_name}_output")
    mv ${current_dir}/clinfo/*cl_output ${current_dir}/cl_sort/
  done

  ### Tree building and log likelihood collection ###
  cd ${current_dir}
  for aln_data in ${current_dir}/Permutation_dataset/*; do
    data_name=$(basename ${aln_data} | cut -d "." -f 1)

    for partition in ${current_dir}/cl_sort/*.cl_output; do
      k_number=$(basename ${partition} | cut -d "." -f 1)

      iqtree -s ${aln_data} -S ${partition} -m ${single_model} --prefix test-temp -c AUTO -seed ${iqtree_seed_number} -keep-ident -quiet 2>&1 >/dev/null

      loglkh=$(cat test-temp.iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
      echo -e "${data_name}\t${k_number}\t${loglkh}" >>${current_dir}/loglkh.collection
      rm -rf ${current_dir}/test-temp.*

    done

  done

  mv ${current_dir}/Permutation_dataset/actual_data.fasta ${current_dir}/aln_concatenation.fasta
  rm -rf ${current_dir}/Permutation_dataset/

  ### Graph plotting for the result ###
  Rscript <(echo "${dloglkhPlot}") 2>&1 >/dev/null
  mv ${current_dir}/permutation_test.svg ${work_dir}/Analysis/

else
  cd ${work_dir}/Analysis/optimal_k/
  current_dir=$(pwd)
  clustree='suppressPackageStartupMessages(library(clustree))
    dd <- read.table("cluster.summary", header = TRUE)
    p <- clustree(dd, prefix = "k"); pdf(NULL)
    ggsave("cluster_summary.svg", width = 12, height = 12)'
  Rscript <(echo "${clustree}") 2>&1 >/dev/null

  ### Visualize the cluster ###
  decision_index=1
  while [[ "${decision_index}" -le 1 ]]; do
    if [[ -e ${work_dir}/Analysis/cluster.visualization ]]; then
      rm -f ${work_dir}/Analysis/cluster.visualization
    fi

    read -p $'\nCheck the cluster_summary.svg and specify the number of cluser to visualize: ' k_visualization
    k_visualization=$(printf "%02d" "${k_visualization}")
    k_visualization_mod=$(echo "k${k_visualization}.cl")

    echo "genes M1 M2 cluster" >${work_dir}/Analysis/cluster.visualization
    while read line; do
      name=$(echo "${line}" | cut -d " " -f 1)
      k=$(cat ${current_dir}/clinfo/${k_visualization_mod} | grep "${name}" | cut -d " " -f 2)
      echo "${line} ${k}" >>${work_dir}/Analysis/cluster.visualization
    done <${work_dir}/Analysis/2-MDS.coord

    cd ${work_dir}/Analysis/
    scatter_plot='suppressPackageStartupMessages(library(ggplot2)) 
      dd <- read.table("cluster.visualization", header = TRUE)
      dd$cluster <- factor(dd$cluster)
      p <- ggplot(data = dd, aes(x = M1,y = M2, color = cluster)) + geom_point(pch = 1, size = 3)
      pdf(NULL)
      ggsave("cluster_visualization.svg", width = 12, height = 12)'
    Rscript <(echo "${scatter_plot}") 2>&1 >/dev/null

    read -p $'\nCheck cluster_visualization.svg.\nAre you happy with the number of cluster?(y/n) ' k_visualization_decision
    if [[ "${k_visualization_decision}" == "y" ]] || [[ "${k_visualization_decision}" == "Y" ]]; then
      ((decision_index++))
    fi
  done
fi
((step++))
echo -e "\n--------------------"

# Part 4. Post analysis operation #
###################################
echo "
Part ${step}. Post analysis operation"
if [[ "${infer_optimal_k}" == "clustree" ]]; then
  cd ${current_dir}
  mkdir ${current_dir}/cl_sort
  potential_k=$(echo "k${k_visualization}")

  PartitionSorting ${current_dir} ${potential_k}.cl
  mv ${current_dir}/clinfo/${potential_k}.cl_output ${current_dir}/cl_sort/
  mv ${work_dir}/aln_concatenation.fasta ${current_dir}
else
  read -p "
  Specify the optimal number of clusters 
  for supergene tree reconstruction: " potential_k
  if [[ -n "${potential_k}" ]] && [[ "${potential_k}" -eq "${potential_k}" ]] 2>/dev/null; then
    potential_k=$(echo "k$(printf "%02d" "${potential_k}")")
  else
    echo -e "\nMessage:The number of clusters failed to be specified."
    exit
  fi
fi

## Supergene tree construction ##
partition_file_name=$(echo "${potential_k}.cl_output")
iqtree -s aln_concatenation.fasta -S ${current_dir}/cl_sort/${partition_file_name} -m ${single_model} -c AUTO -seed ${iqtree_seed_number} -keep-ident -quiet 2>&1 >/dev/null
mv ${current_dir}/cl_sort/${partition_file_name}.parstree ${current_dir}/clusters_topology.nw
mv ${current_dir}/cl_sort/${partition_file_name}.treefile ${current_dir}/clusters_phylo.nw
rm -rf ${current_dir}/clinfo/ ${current_dir}/cl_sort/ ${current_dir}/partition
mv ${current_dir}/* ${work_dir}/Analysis/
rm -rf ${current_dir}

((step++))
time_stamp=$(date)
echo -e "\n--------------------\nAnalysis finished at\n${time_stamp}."
