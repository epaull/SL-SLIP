#!/usr/bin/env Rscript

cc <- as.numeric(as.matrix(read.table("similarity.LIN-RESNIK.CC.subset.tab"))[,3])
mf <- as.numeric(as.matrix(read.table("similarity.LIN-RESNIK.MF.subset.tab"))[,3])
bp <- as.numeric(as.matrix(read.table("similarity.LIN-RESNIK.BP.subset.tab"))[,3])

png("mf-before.png")
plot(density(mf),ylim=c(0,1))
png("mf-after.png")
plot(density(mf^0.33),ylim=c(0,1))


png("cc-before.png")
plot(density(cc),ylim=c(0,1))
png("cc-after.png")
plot(density(cc^0.33),ylim=c(0,1))

png("bp-before.png")
plot(density(bp),ylim=c(0,1))
png("bp-after.png")
plot(density(bp^0.33),ylim=c(0,1))

q();
