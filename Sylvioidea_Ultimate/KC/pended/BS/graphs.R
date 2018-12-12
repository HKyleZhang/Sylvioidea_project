# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)

out.file <- ""
file.names <- dir(pattern = "logLkh")
for(i in 1:length(file.names)){
  file<-read.table(file.names[i])
  out.file <- cbind(out.file, file)
}

data <- out.file[,2:16]
data <- t(data)
write.csv(data, file = "logLkhcollection.csv", row.names = FALSE)
data <- read.csv("logLkhcollection.csv")
k <- seq(1:15)
Mean <- apply(data, 1, mean)
SD <- apply(data, 1, sd)

DeltaK <- data.frame(Mean)
for (i in 1:14){
  DeltaK[i,] <- DeltaK[i,] - DeltaK[i+1,]
}

data <- cbind(k, DeltaK, Mean, SD, data)
num <- seq(1:100)
names(data) <- c("K", "Delta_K", "Mean", "SD", num)
data$Mean <- 0-data$Mean

p1 <- ggplot(data,aes(x=K,y=Mean)) + geom_point() + geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD))
p2 <- ggplot(data[1:14,],aes(x=data[1:14,1],y=data[1:14,2])) + geom_line() + geom_point()
