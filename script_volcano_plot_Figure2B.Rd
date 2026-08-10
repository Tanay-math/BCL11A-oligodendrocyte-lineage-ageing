
## volcano plot

setwd("path/to/directory")

D<-read.table("GSE303317_Differential_expression_agedOPC_vs_youngOPC.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)  #  GSE303317_Differential_expression_agedOPC_vs_youngOPC.txt  is available at GSE303317

D_na_omit<-na.omit(D)
up<-D_na_omit[,c(8)]=="UP"
D_na_omit_up<-D_na_omit[up,]
print("dim of up-regulated list")
print(dim(D_na_omit_up))
down<-D_na_omit[,c(8)]=="DOWN"
D_na_omit_down<-D_na_omit[down,]
print("dim of down-regulated list")
print(dim(D_na_omit_down))

par(mar=c(5.1,5.1,4.1,2.1))

plot(D_na_omit$log2FoldChange,-log(D_na_omit$padj,10) ,type="p", pch=21, bg="#96969655", cex=0.35, lwd=0,  xlim=c(-10,10), ylim=c(0,8), bty="n", main="", xlab=expression(paste(log[2],"(Aged/Young)")), ylab=expression(paste(-log[10], "(padj)")))

abline(v=c(-0.58,0.58), col="magenta",  lty=2, lwd=1)
abline(h=1.33, col="blue",  lty=2, lwd=1)
par(new=TRUE)
plot(D_na_omit_up$log2FoldChange,-log(D_na_omit_up$padj,10) ,type="p", pch=21, bg="red", cex=0.35, lwd=0,xlim=c(-10,10), ylim=c(0,8), axes=FALSE, main="", xlab="", ylab="", bty="n")
par(new=TRUE)
plot(D_na_omit_down$log2FoldChange,-log(D_na_omit_down$padj,10) ,type="p", pch=21, bg="green", cex=0.35, lwd=0,xlim=c(-10,10), ylim=c(0,8), axes=FALSE, main="", xlab="", ylab="", bty="n")

print("success volcano plot")




