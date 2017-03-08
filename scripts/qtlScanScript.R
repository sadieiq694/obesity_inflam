library(qtl)

load("data/BTBR.clean.data.Rdata")
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
my.phenos <- phenotypes.rz[,c("WT.10wk", "Fat.wt", "Weight")]

my.phenos
head(my.phenos)
f2g$pheno <- cbind(f2g$pheno, my.phenos)
names(f2g$pheno)

f2g<- calc.genoprob(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)
f2g <- sim.geno(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)

sex <- as.numeric(f2g$pheno$Sex)
f2g.scan1 <-scanone(f2g,  pheno.col = 4:6, addcovar = sex, method = "hk")

# f2g.perm1 <-scanone(f2g, pheno.col = 4:6, addcovar = sex, method= "hk", n.perm=100, perm.Xsp = TRUE)
# save(list="f2g.perm1", file="f2g_perm1.Rdata")

load(file="data/f2g_perm1.Rdata")

summary(f2g.scan1, threshold=3)
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#my.first.scan <- scanone(f2g, pheno.col = 4:6, addcovar = sex,
#                         method = "hk") #b/c first three are identifiers etc

#my.first.perm <- scanone(f2g, pheno.col = 4:6, addcovar = sex,
#                         method = "hk", n.perm = 100, perm.Xsp = TRUE)
for(i in 1:3){
  plot(f2g.scan1, lodcolumn = i)
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
}

# add gene as a covariate for a qtl scan??

lodint(f2g.scan1, chr = 10)

?calc.genoprob()
