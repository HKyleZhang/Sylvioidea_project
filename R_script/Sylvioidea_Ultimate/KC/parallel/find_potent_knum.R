logd <- read.csv("logLkhforR.csv", header = TRUE)
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
write.table(res, file = "potent_knum", quote = FALSE, col.names = FALSE, row.names = FALSE)
