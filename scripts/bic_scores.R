#USE FOR SHARED QTL PEAKS?

rm(list=ls())
setwd("C:/Users/sadie.la/Documents/obesity_inflam")

load("data/BTBR.clean.data.Rdata")

phenotypes.rz$Fat.wt[phenotypes.rz$Fat.wt <0 & is.numeric(phenotypes.rz$Fat.wt)] <- NA

library(qtl)
library(ggplot2)
library(qtlnet)

load("data/BTBR.clean.data.Rdata")
names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]

Il1b.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.islet<- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

Il1b.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Il1b")]]
Nfkb1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nfkb1")]]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("Fat.wt", "Weight", "adipose.turnover")], Il1b.islet, Nfkb1.islet, Il1b.adipose, Nfkb1.adipose)
names(f2g$pheno)

#source("important_func.R")

#####################################
#       Functions                   #
#####################################

#get gene expression data for a gene with a given name from a tissue
#gene name: MGI gene symbol
#data.set: gene expression dataset, such as adipose.rz
gene.exp <- function(gene.name, data.set) {
  return(data.set[,annot$a_gene_id[which(annot$gene_symbol==gene.name)]])
}


#Get clinical data for a certain parameter
#clin.name: name of the clinical trait in the dataset
#data.set: either  phenotypes or phenotypes.rz
clinical <- function(clin.name, data.set) {
  if(clin.name %in% names(data.set)) {
    return(data.set[,clin.name])
  }
}

#Get the genotype of the SNPs at a certain position on a chromosome
#chr: the chromosome on which to look for the genotype data
#pos: the position at which to look for genotyp data
genotype <- function(chr, pos) {
  if(chr >0 && chr < 21 && pos >= 0 && pos <= 100) {
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr = chr, pos= pos)])
  }
}

#BIC Score model model analysis
#X: gene expression data for a given gene
#Y: quantitative measurement of clinical phenotype
#Q: genotype at a given marker
triple.fit <- function(X, Y, Q) {
  
  #Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(X)), which(is.na(Y)), which(is.na(Q)))))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  print(paste("Removed", length(indx), " rows with NA values from data.", sep = ""))

  #Calculate BIC scores for models
  bic.independent <- BIC(lm(X~Q)) + BIC(lm(Y~Q)) #X<-Q->Y
  bic.reactive <- BIC(lm(X~Y)) + BIC(lm(Y~Q)) # Q->Y->X
  bic.causal <- BIC(lm(X~Q)) + BIC(lm(Y~X)) # Q->X->Y
  bic.complex <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))
  
  # Print out the scores from each model
  print("BIC Scores of each model")
  scores <- c(bic.independent, bic.reactive, bic.causal, bic.complex)
  names(scores) <- c("independent", "reactive", "causal", "complex")
  print(scores)
  
  # Make lowest BIC score 0 and linearize all other scores accordingly to calculate Delta values
  deltas <- scores - min(scores)
  
  # Estimate the strength of evidence for each model
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  # Print out the probabilities of each model being the likely explanation for the data
  print("Probability of each model explaining the data")
  print(strengths * 100)
  
  # Print out how many more times likely the best model is
  print("The factor by which the best model is better than the rest")
  print(max(strengths) / strengths)
  
}

#peak on chr 2 at position 70
triple.fit(gene.exp("Il1b", islet.rz), clinical("Fat.wt", phenotypes.rz), genotype(chr = 2, pos = 70))
triple.fit(gene.exp("Il1b", islet.rz), clinical("Fat.wt", phenotypes.rz), genotype(chr = 2, pos = 73.7))

C12 <- genotype(chr = 12, pos = 24.4)


triple.fit(Nfkb1.islet, phenotypes.rz$adipose.turnover, C12)

triple.fit(gene.exp("Nfkb1", islet.rz), clinical("Adipose.turnover", phenotypes.rz), genotype(chr = 12, pos = 24.4))
#peak chromosome 10 pos 42
#peak chromosome 12 pos 24.4

# fat weight: chrom 2, pos 70
# body weight: chrom 2 pos 56.3
# Il1b islet expression: chrom 2 pos 73.7








