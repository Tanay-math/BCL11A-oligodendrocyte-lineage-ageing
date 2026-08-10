#making PCA
setwd("path/to/directory")
library(PCAtools)
dat <- read.table("normalised_lineage_data.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)  #  normalised_lineage_data is available at GSE301578
dim(dat)
#because of 0 values, we do dat+1
dat_log<-log(dat+1,2)
a<-rowSums(dat_log[,c(1:12)]==0)==12
dat_log_a<-dat_log[a,]
dim(dat_log_a)
b<-setdiff(rownames(dat_log),rownames(dat_log_a))
dat_log_b<-dat_log[b,]
dim(dat_log_b)
metadata <- data.frame(row.names = colnames(dat_log_b))
metadata$Group <- rep(NA, ncol(dat_log_b))
metadata$Group[seq(1,4,1)]<-'OPC'
metadata$Group[seq(5,8,1)]<-'PreOL'
metadata$Group[seq(9,12,1)]<-'OL'
p<-pca(dat_log_b,scale=T, center=T, metadata=metadata)
biplot(p, colby='Group', colkey=c(OPC="red",PreOL="green", OL="blue"),legendPosition = "right",hline=-30,vline=-25, title='PCA bi-plot', drawConnectors=FALSE, lab = NULL, gridlines.major = FALSE, gridlines.minor = FALSE)

#making cluster
library(WGCNA)
sampleTree = hclust(dist(t(dat_log_b)), method = "average");
plot(sampleTree, main = "Sample clustering ", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

