library(gplots)
setwd("path/to/Data")
##type3
Dat<-read.table("output_type3_all_affected_in_age.csv", sep=",", header=T, row.names=1)
D<-Dat[-c(1:7)]

heatmap.2(as.matrix(D), col=greenred(75), scale="row", key=T, keysize=1.5, density.info="none", trace="none",cexCol=0.9,Colv=FALSE)


##type4
Dat<-read.table("output_type4_all_genes_affected_in_age.csv", sep=",", header=T, row.names=1)
D<-Dat[-c(1:7)]

heatmap.2(as.matrix(D), col=greenred(75), scale="row", key=T, keysize=1.5, density.info="none", trace="none",cexCol=0.9,Colv=FALSE)


## heatmap_5TFs

Dat<-read.table("for_heatmap_5TFs.txt", sep="\t", header=T, row.names=1)

heatmap.2(as.matrix(Dat), col=greenred(75), scale="row", key=T, keysize=1.5, density.info="none", trace="none",cexCol=0.9,Colv=FALSE)



