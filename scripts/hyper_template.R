rm(list=ls())

setwd("C/Users/sadie.la/Documents/obesity_inflam")
library(qtl)

load("data/BTBR.clean.data.Rdata")
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]

Il1b.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

Il1b.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("Fat.wt", "Weight", "adipose.turnover")], Il1b.islet, Nfkb1.islet, Il1b.adipose, Nfkb1.adipose)
names(f2g$pheno)

f2g<- calc.genoprob(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)
f2g <- sim.geno(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)

sex <- as.numeric(f2g$pheno$Sex)

load(file="data/f2g_perm1.Rdata")

genotype <- function(chr, pos) {
  if(chr >0 && chr < 21 && pos >= 0 && pos <= 100) {
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr = chr, pos= pos)])
  }
}


C2 <- genotype(2, 70)
C2
  
# fat weight: chrom 2, pos 70
# body weight: chrom 2 pos 56.3
# Il1b islet expression: chrom 2 pos 73.7

phenotypes.rz$Fat.wt[phenotypes.rz$Fat.wt <0 & is.numeric(phenotypes.rz$Fat.wt)] <- NA


#Shared QTL Peak Chrom 2 Fat Weight, Il1b islet expression
#LOW SCORE IS BEST!
BIC(lm(formula = Fat.wt ~ sex, data = f2g$pheno))
#1392.893
BIC(lm(formula = Fat.wt ~ sex + Il1b.islet, data = f2g$pheno))
#1268.208
BIC(lm(formula = Fat.wt ~ sex + Il1b.islet + C2, data = f2g$pheno))
#1239.495




#Shared QTL Peak chrom 12 adipose turnover, Il1b



f2g$pheno <- transform(f2g$pheno, Q2 = as.factor(f2g$geno[[]]$data[,find.marker(f2g, CHR, POS)]))
levels(f2g$pheno$Q2) <- c("B", "H", "R")
#f2g <- fill.geno(f2g, method="argmax")
#f2g$pheno <- transform(f2g$pheno,Q2 = as.factor(f2g$geno[[2]]$data[,find.marker(f2g, 2, 75.2)]))
#levels(f2g$pheno$Q2) <- c("B","H","R")names(f2g$pheno)f2g$pheno <- transform(f2g$pheno, adipor1_islet)
qplot(adipor1_islet, INS.10wk, color=Q2, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q2), method="lm", se=FALSE)
#Ghrh_islet <- islet.rz[, annot[grep("Ghrh$", annot$gene1), 1]]
#lm(formula = INS.10wk ~ Adipor1_islet + Q2, data = f2g$pheno)
# ^ do summary of this