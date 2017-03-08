library(qtl)

load("data/BTBR.clean.data.Rdata")
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]

Il1b.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]
Il1b.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

#THIS LINE DOESNT WORK
#my.phenos <- phenotypes.rz[,c("WT.10wk", "Fat.wt", "Weight", Il1b.islet, Nfkb1.islet, Il1b.kidney, Nfkb1.kidney, Il1b.adipose, Nfkb1.adipose, Il1b.liver, Nfkb1.liver)]


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("Fat.wt", "Weight", "adipose.turnover")],Il1b.islet, Nfkb1.islet, Il1b.adipose, Nfkb1.adipose)


#my.phenos
#head(my.phenos)
# f2g$pheno <- cbind(f2g$pheno, my.phenos)
names(f2g$pheno)

f2g<- calc.genoprob(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)
f2g <- sim.geno(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)

sex <- as.numeric(f2g$pheno$Sex)
f2g.scan1 <-scanone(f2g,  pheno.col = 4:10, addcovar = sex, method = "hk")

# f2g.perm1 <-scanone(f2g, pheno.col = 4:6, addcovar = sex, method= "hk", n.perm=100, perm.Xsp = TRUE)
# save(list="f2g.perm1", file="f2g_perm1.Rdata")

load(file="data/f2g_perm1.Rdata")

summary(f2g.scan1, threshold=3)

for(i in 1:7){
  plot(f2g.scan1, lodcolumn = i)
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
}

# add gene as a covariate for a qtl scan??

lodint(f2g.scan1, chr = 10)

