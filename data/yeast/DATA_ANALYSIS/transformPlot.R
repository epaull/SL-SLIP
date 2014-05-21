


mf <- as.numeric(as.matrix(read.table("similarity.LIN-RESNIK.MF.subset.tab"))[,3])
bp <- as.numeric(as.matrix(read.table("similarity.LIN-RESNIK.BP.subset.tab"))[,3])
plot(density(mf^0.33),ylim=c(0,1))
plot(density(cc^0.33),ylim=c(0,1))
plot(density(bp^0.33),ylim=c(0,1))

