#Differential expression analysis OPC young vs OPC aged

#Batcheffect removed from the data and PCA to see.

setwd("path/to/directory")
library(tximport);
library(DESeq2);
print("library ready")
print("reading files")

files<-c("GSM9085542_OPC-quant.sf","GSM9085543_OPC-quant.sf","GSM9085544_OPC-quant.sf", "GSM9085545_OPC-quant.sf", "GSM3967160_OPC-quant.sf", "GSM3967161_OPC-quant.sf", "GSM3967162_OPC-quant.sf", "GSM3967163_aged-quant.sf", "GSM3967164_aged-quant.sf","GSM3967165_aged-quant.sf", "GSM3929171_aged-quant.sf", "GSM3929172_aged-quant.sf","GSM3929170_aged-quant.sf")

# quant.sf files are available in the Data folder.

names(files)<-c("OPC.1","OPC.2","OPC.3","OPC.4","young1", "young2", "young3","Aged.1","Aged.2","Aged.3","Aged.5","Aged.6","Aged.4")


tx2gene<-read.table("Rnor_mart_export_tx2gene.txt", sep="\t")
txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")


group<-factor(c("OPC","OPC","OPC","OPC","OPC", "OPC", "OPC", "AgeOPC","AgeOPC","AgeOPC","AgeOPC","AgeOPC", "AgeOPC"))
df<-data.frame(group);
rownames(df)<-colnames(txi$counts);

df$batch<-factor(c("I","I","I","I","II","II","II","II","II","II","III","III", "III"))

dds<-DESeqDataSetFromTximport(txi, df, design= ~batch+group)


dds1 = DESeq(dds)

res_effectOfAgeing<-results(dds1, contrast=c("group","AgeOPC","OPC"))

id2s<-read.table("ID2Symbol.txt", sep="\t", header=T, row.names=1)
res_effectOfAgeing_id2s<-as.data.frame(res_effectOfAgeing)
a<-intersect(row.names(res_effectOfAgeing_id2s),row.names(id2s))
id2s_a<-as.data.frame(id2s[a,])
names(id2s_a)<-"Gene_symbol"
res_effectOfAgeing_id2s<-data.frame(id2s_a, res_effectOfAgeing_id2s)



print("Differential expression Successfully completed")

# to see the if batch effect removed
print("removing batch effect")
vsd<-varianceStabilizingTransformation(dds1)
mm <- model.matrix(~group, colData(vsd))
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=df$batch, design=mm)
dat_vsd<-data.frame(id2s_a, assay(vsd))


# PCA to check batch effect
library(PCAtools)

metadata <- data.frame(row.names = colnames(dat_vsd[-1]))
metadata$Group <- rep(NA, ncol(dat_vsd[-1]))
metadata$Group[seq(1,7,1)]<-'Young'
metadata$Group[seq(8,13,1)]<-'Aged'
p<-pca(dat_vsd[-1],scale=T, center=T, metadata=metadata)

biplot(p, colby= 'Group', colkey=c(Young="blue",Aged="black"),legendPosition = "right", hline=-10, title='PCA bi-plot', drawConnectors=FALSE,  lab = NULL, gridlines.major = FALSE, gridlines.minor = FALSE)

print("success")
