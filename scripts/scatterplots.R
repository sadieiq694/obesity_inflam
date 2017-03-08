library(qtl)
library(ggplot2)

load("data/BTBR.clean.data.Rdata")
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]

summary(f2g)
names(f2g)
class(f2g$geno)
names(f2g$geno)
class(f2g$pheno)
names(f2g$pheno)


Il1b.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]
Il1b.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("Fat.wt", "Weight", "adipose.turnover")],Il1b.islet, Nfkb1.islet, Il1b.adipose, Nfkb1.adipose)
names(f2g$pheno)


#scatter plots#

#phenotype distributions
qplot(Fat.wt, facets = Sex~., data = f2g$pheno)
qplot(Weight, facets = Sex~., data = f2g$pheno)
qplot(adipose.turnover, facets = Sex~., data = f2g$pheno)

#gene expression distributions
qplot(Il1b.islet, facets = Sex~., data = f2g$pheno)
qplot(Nfkb1.islet, facets = Sex~., data = f2g$pheno)
qplot(Il1b.adipose, facets = Sex~., data = f2g$pheno)
qplot(Nfkb1.adipose, facets = Sex~., data = f2g$pheno)

#phenotypes against each other
qplot(Fat.wt, Weight, color=Sex, data = f2g$pheno) + geom_smooth(method = "lm") + aes(x= Fat.wt, y=Weight, color=Sex)
qplot(Fat.wt, adipose.turnover, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(Weight, adipose.turnover, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")

#Fat.wt vs gene expressions
#R = -0.25
qplot(Fat.wt, Il1b.islet, color=Sex, data = f2g$pheno) + geom_smooth(method = "lm")
#R = 0.26
qplot(Fat.wt, Nfkb1.islet, color=Sex, data = f2g$pheno) + geom_smooth(method = "lm")
#R = 0.24
qplot(Fat.wt, Il1b.adipose, color=Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(Fat.wt, Nfkb1.adipose, color=Sex, data = f2g$pheno) + geom_smooth(method = "lm")

#Body weight vs. gene expressions
#R = -0.34
qplot(Weight, Il1b.islet, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
#R = 0.36
qplot(Weight, Nfkb1.islet, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(Weight, Il1b.adipose, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(Weight, Nfkb1.adipose, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")

#adipose turnover vs. gene expressions
qplot(adipose.turnover, Il1b.islet, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
#R = 0.36
qplot(adipose.turnover, Nfkb1.islet, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(adipose.turnover, Il1b.adipose, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(adipose.turnover, Nfkb1.adipose, color = Sex, data = f2g$pheno)

#gene expressions from same tissue
qplot(Nfkb1.islet, Il1b.islet, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")
qplot(Nfkb1.adipose, Il1b.adipose, color = Sex, data = f2g$pheno) + geom_smooth(method = "lm")




