#!/usr/bin/env Rscript

aucs <- as.matrix(read.table("aucs.txt",sep="\t",header=TRUE))[,2:9]
pdf("aucs.pdf")
plot(aucs[,3],aucs[,2],type='p',col='blue',lwd=2,main="AUC for collective inference (Red) vs non-collective (Blue)")
lines(aucs[,3],aucs[,1],type='p',col='red',lwd=2)
abline(h=0.5, lty=2)

pdf("params.pdf")
plot(aucs[,9], aucs[,2], ylab='AUC', xlab=colnames(aucs)[9], col='red', type='p',xlim=c(0,10))


q();
