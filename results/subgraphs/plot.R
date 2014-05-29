#!/usr/bin/env Rscript

aucs <- as.matrix(read.table("aucs.txt"))[,2:4]
pdf("aucs.pdf")
plot(aucs[,3],aucs[,2],type='p',col='blue',lwd=2,main="AUC for collective inference (Red) vs non-collective (Blue)")
lines(aucs[,3],aucs[,1],type='p',col='red',lwd=2)
abline(h=0.5, lty=2)
q();
