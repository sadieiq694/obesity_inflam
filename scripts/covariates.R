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

f2g$pheno <- transform(f2g$pheno, Q2 = as.factor(f2g$geno[[]]$data[,find.marker(f2g, CHR, POS)]))
levels(f2g$pheno$Q2) <- c("B", "H", "R")
#f2g <- fill.geno(f2g, method="argmax")
#f2g$pheno <- transform(f2g$pheno,Q2 = as.factor(f2g$geno[[2]]$data[,find.marker(f2g, 2, 75.2)]))
#levels(f2g$pheno$Q2) <- c("B","H","R")names(f2g$pheno)f2g$pheno <- transform(f2g$pheno, adipor1_islet)
qplot(adipor1_islet, INS.10wk, color=Q2, shape=Sex, data=f2g$pheno) + geom_smooth(aes(group=Q2), method="lm", se=FALSE)
#Ghrh_islet <- islet.rz[, annot[grep("Ghrh$", annot$gene1), 1]]
#lm(formula = INS.10wk ~ Adipor1_islet + Q2, data = f2g$pheno)
# ^ do summary of this

#SHARED QTL PEAKS:
#fat weight and Il1b islet expression, chromosome 2
#body weight and Il1b islet expression, chromosome 2
#adipose turnover and Il1b islet expression, chromosome 12
#body weight Nfkb1 islet 5

#fat weight
f2g.scanFat <- scanone(f2g, pheno.col = 4, addcovar = sex, method= "hk")
max(f2g.scanFat)
plot(f2g.scanFat)
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanFat, chr = 2) #<-- can specify chromosome for which you want a peak! (use the same chr argument for plots)
max(f2g.scanFat, chr = 10)
#peak on chr 2 at position 70
#peak on chr 10 at position 48.3


#fat weight with Il1b islet expression covariate
#significant qtl drop chromosome 2
f2g.Fat.Il1b.islet.covar <- scanone(f2g, pheno.col = 4, addcovar = Il1b.islet, method= "hk")
plot(f2g.Fat.Il1b.islet.covar)
plot(f2g.scanFat, f2g.Fat.Il1b.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#fat weight with Nfkb1 islet expression covariate
f2g.Fat.Nfkb1.islet.covar <- scanone(f2g, pheno.col = 4, addcovar = Nfkb1.islet, method= "hk")
plot(f2g.Fat.Nfkb1.islet.covar)
plot(f2g.scanFat, f2g.Fat.Nfkb1.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#fat weight with Il1b adipose expression covariate
f2g.Fat.Il1b.adip.covar <- scanone(f2g, pheno.col = 4, addcovar = Il1b.adipose, method= "hk")
plot(f2g.Fat.Il1b.adip.covar)
plot(f2g.scanFat, f2g.Fat.Il1b.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#fat weight with Nfkb1 adipose expression covariate
f2g.Fat.Nfkb1.adip.covar <- scanone(f2g, pheno.col = 4, addcovar = Nfkb1.adipose, method= "hk")
plot(f2g.Fat.Nfkb1.adip.covar)
plot(f2g.scanFat, f2g.Fat.Nfkb1.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanFat, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Body weight
f2g.scanWeight <- scanone(f2g, pheno.col = 5, addcovar = sex, method= "hk")
max(f2g.scanWeight)
#max(f2g.scanWeight, chr = c(1,3,15)) <-- can specify chromosome for which you want a peak! (use the same chr argument for plots)
plot(f2g.scanWeight)
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanWeight, chr = 2)
max(f2g.scanWeight, chr = 5)
max(f2g.scanWeight, chr = 16)
#peak chromosome 2 position 56.3
#peak chromosome 5 position 54.7
#peak chromosome 16 position 39.2

#Body weight with Il1b islet expression covariate
f2g.Weight.Il1b.islet.covar <- scanone(f2g, pheno.col = 5, addcovar = Il1b.islet, method= "hk")
plot(f2g.Weight.Il1b.islet.covar)
plot(f2g.scanWeight, f2g.Weight.Il1b.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Body weight with Nfkb1 islet expression covariate
#QTL drop on chromosome 5 (pos 54.7)
f2g.Weight.Nfkb1.islet.covar <- scanone(f2g, pheno.col = 5, addcovar = Nfkb1.islet, method= "hk")
plot(f2g.Weight.Nfkb1.islet.covar)
plot(f2g.scanWeight, f2g.Weight.Nfkb1.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Body weight with Il1b adipose expression covariate
f2g.Weight.Il1b.adip.covar <- scanone(f2g, pheno.col = 5, addcovar = Il1b.adipose, method= "hk")
plot(f2g.Weight.Il1b.adip.covar)
plot(f2g.scanWeight, f2g.Weight.Il1b.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Body weight with Nfkb1 adipose expression covariate
f2g.Weight.Nfkb1.adip.covar <- scanone(f2g, pheno.col = 5, addcovar = Nfkb1.adipose, method= "hk")
plot(f2g.Weight.Nfkb1.adip.covar)
plot(f2g.scanWeight, f2g.Weight.Nfkb1.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanWeight, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Adipose Turnover
f2g.scanAdTurn <- scanone(f2g, pheno.col = 6, addcovar = sex, method= "hk")
max(f2g.scanAdTurn)
plot(f2g.scanAdTurn)
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanAdTurn, chr= 1)
max(f2g.scanAdTurn, chr= 10)
max(f2g.scanAdTurn, chr= 12)
#peak chromosome 1 pos 89
#peak chromosome 10 pos 42
#peak chromosome 12 pos 24.4

#Adipose Turnover with Il1b islet expression covariate
#significant QTL drop chromosome 10, 12
f2g.Ad.Turn.Il1b.islet.covar <- scanone(f2g, pheno.col = 6, addcovar = Il1b.islet, method= "hk")
plot(f2g.Ad.Turn.Il1b.islet.covar)
plot(f2g.scanAdTurn, f2g.Ad.Turn.Il1b.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Adipose Turnover with Nfkb1 islet expression covariate
#significant QTL drop chromosome 10
f2g.Ad.Turn.Nfkb1.islet.covar <- scanone(f2g, pheno.col = 6, addcovar = Nfkb1.islet, method= "hk")
plot(f2g.Ad.Turn.Nfkb1.islet.covar)
plot(f2g.scanAdTurn, f2g.Ad.Turn.Nfkb1.islet.covar, col=c("black", "red"))
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Adipose Turnover with Il1b adipose expression covariate
#significant QTL drop chromosomes 10, 12
f2g.Ad.Turn.Il1b.adip.covar <- scanone(f2g, pheno.col = 6, addcovar = Il1b.adipose, method= "hk")
plot(f2g.Ad.Turn.Il1b.adip.covar)
plot(f2g.scanAdTurn, f2g.Ad.Turn.Il1b.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Adipose Tissue with Nfkb1 adipose expression covariate
#significant QTL drop chromosomes 10, 12
#MOST SIGNIFICANT ADIPOSE TURNOVER DROP
f2g.Ad.Turn.Nfkb1.adip.covar <- scanone(f2g, pheno.col = 6, addcovar = Nfkb1.adipose, method= "hk")
plot(f2g.Ad.Turn.Nfkb1.adip.covar)
plot(f2g.scanAdTurn, f2g.Ad.Turn.Nfkb1.adip.covar, col=c("black", "red"))
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanAdTurn, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Il1b islet expression
f2g.scanIlIs <- scanone(f2g, pheno.col = 7, addcovar = sex, method= "hk")
max(f2g.scanIlIs)
plot(f2g.scanIlIs)
add.threshold(f2g.scanIlIs, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanIlIs, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanIlIs, chr = 1)
max(f2g.scanIlIs, chr = 2)
max(f2g.scanIlIs, chr = 12)
max(f2g.scanIlIs, chr = 13)
#peak chromosome 1 pos 31.4
#peak chromosome 2 pos 73.7
#peak chromosome 12 pos 29.8
#peak chromosome 13 pos 68.7

#Il1b adipose expression
#located on chromosome 2 location 129364570-129371139
f2g.scanIlAd <- scanone(f2g, pheno.col = 9, addcovar = sex, method= "hk")
max(f2g.scanIlAd)
plot(f2g.scanIlAd)
add.threshold(f2g.scanIlAd, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanIlAd, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanIlAd, chr = 13)
max(f2g.scanIlAd, chr = 14)
#peak chromosome 13 pos 63
#peak chromosome 14 pos 63.9

#Checking for shared qtls: fat weight and Il1b gene expression
#shared peak chromosome 2 pos 73.7
plot(f2g.scanFat, f2g.scanIlIs, f2g.scanIlAd, col=c("black", "red", "blue"))

#Checking for shared qtls: body weight and Il1b gene expression
#shared peak chromosome 2?
plot(f2g.scanWeight, f2g.scanIlIs, f2g.scanIlAd, col=c("black", "red", "blue"))

#Checking for shared qtls: adipose turnover and Il1b gene expression
plot(f2g.scanAdTurn, f2g.scanIlIs, f2g.scanIlAd, col=c("black", "red", "blue"))

#Nfkb1 islet expression
f2g.scanNfIs <- scanone(f2g, pheno.col = 8, addcovar = sex, method= "hk")
max(f2g.scanNfIs)
plot(f2g.scanNfIs)
add.threshold(f2g.scanNfIs, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanNfIs, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanNfIs, chr = 3)
max(f2g.scanNfIs, chr = 6)
max(f2g.scanNfIs, chr = 17)
#peak chromosome 3 pos 63
#peak chromosome 6 pos 91.4
#peak chromosome 17 pos 19

#Nfkb1 adipose expression
f2g.scanNfAd <- scanone(f2g, pheno.col = 10, addcovar = sex, method= "hk")
max(f2g.scanNfAd)
plot(f2g.scanNfAd)
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
summary(f2g.scanFat, threshold=3)
max(f2g.scanNfAd, chr = 6)
max(f2g.scanNfAd, chr = 17)
#peak chromosome 6 pos 67.9
#peak chromosome 17 pos 11.8

#Checking for shared qtls: fat weight and Nfkb1 gene expression
plot(f2g.scanFat, f2g.scanNfIs, f2g.scanNfAd, col=c("black", "red", "blue"))
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
#shared peak fat, Nfkb1 islet expression, chr 5?

#Checking for shared qtls: body weight and Nfkb1 gene expression
plot(f2g.scanWeight, f2g.scanNfIs, f2g.scanNfAd, col=c("black", "red", "blue"))
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

#Checking for shared qtls: adipose turnover and Nfkb1 gene expression
plot(f2g.scanAdTurn, f2g.scanNfIs, f2g.scanNfAd, col=c("black", "red", "blue"))
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scanNfAd, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

