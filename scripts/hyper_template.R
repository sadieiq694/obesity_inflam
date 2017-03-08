rm(list=ls())
data(hyper)
hyper <- sim.geno(hyper, step = 2.5, n.draws=16, err=0.01)
out1 <- scanone(hyper, method = "imp")
plot(out1)
max(out1)
find.marker(hyper, 4, 29.5) #marker is "D$Mit164"

#getting rid of peak to better see other QTLs
g <- pull.geno(hyper)[,"D4Mit164"]
out.c4 <- scanone(hyper, method="imp", addcovar=g)
plot(out1, out.c4, col=c("blue", "red"))

#looking for loci that interact with the chr4 locus
out.c4i <- scanone(hyper, method="imp", addcovar=g, intcovar=g)
plot(out.c4i - out.c4)