# Pheno/Gene Expression Reevaluation
# Sadie Allen
# July 11, 2017
# Correlations of more current pheno/gene list

rm(list = ls())

library(corrplot)

# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

names(ghrelin_list)

ghrelin_corr <- ghrelin_list[,5:15]

pcor <- cor(ghrelin_corr, use= "complete.obs")
round(pcor, digits = 2)
corrplot(pcor, order = "hclust")

# efnb3 and ghsr have very similar correlations patterns, but none of the chromosome 18
# candidate genes are strongly correlated with ghsr or efnb3! 
