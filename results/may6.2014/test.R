

vals <- as.matrix(read.table("corr.in"))
ks.test(vals[which(vals[,2] < 0.1),1], vals[which(vals[,2] > 0.1),1])
