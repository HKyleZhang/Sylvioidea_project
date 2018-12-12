library(treespace)

files <- list.files(pattern = "*tre")
for (i in files){
  name <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  name <- paste(name, ".kcd", sep = "")
  trees <- read.tree(i)
  res <- multiDist(trees, lambda = 0.5)
  res <- as.matrix(res)
  kcd <- round(res[2:length(res[,1]),1],8)
  write.table(kcd, file = name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
