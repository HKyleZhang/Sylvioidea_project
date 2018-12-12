args <- commandArgs(TRUE)
library(treespace)
library(foreach)
library(doParallel)

no_cores <- detectCores()
registerDoParallel(no_cores)

L <- as.numeric(args[1])

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

stopImplicitCluster()