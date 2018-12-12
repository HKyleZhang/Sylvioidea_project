library(ggplot2)
library(dplyr)
library(svglite)
data <- read.table("Cluster_Assignment.summary", header = FALSE)
datash <- read.table("Cluster_Assignment_shortdata.summary", header = FALSE)

names(data) <- c("order", "gene","cluster","percentage")
data$percentage <- data$percentage / 10
data$order <- factor(data$order)

names(datash) <- c("order", "gene","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")

res<-ggplot() + geom_bar(data=data, aes(y = percentage,x = order,fill = cluster), stat="identity", position = "stack", width = 0.9) + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, angle = 90), axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0),"cm")) + scale_x_discrete(position = "bottom", labels = datash$gene) + scale_fill_manual(values = c("coral2", "aquamarine4", "yellow3", "deepskyblue2", "magenta2"))

ggsave("Cluster_Assignment.svg", width = 15, height = 5) 
