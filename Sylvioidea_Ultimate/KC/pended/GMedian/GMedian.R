args <- commandArgs(TRUE)

library(treespace)

trees <- read.tree(args[1])

mediantre <- medTree(trees, lambda = 0.5)$treenumbers[1]
mediantre <- as.data.frame(mediantre)
name <- rownames(mediantre)
write.table(name,file = "Median.tre", quote = FALSE, row.names = FALSE, col.names = FALSE)
