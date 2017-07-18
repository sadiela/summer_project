# Pheno/Gene Expression Reevaluation
# Sadie Allen
# July 11, 2017
# Correlations of more current pheno/gene list

rm(list = ls())

library(corrplot)

source("scripts/functions.R")

# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

names(ghrelin_list)

ghrelin_corr <- ghrelin_list[,c(5,6,7,8,9,15)]
arg1_exp <- get_exp_dat("Arg1")
ptprz1_exp <- get_exp_dat("Ptprz1")
kcnd2_exp <- get_exp_dat("Kcnd2")
cacna1h_exp <- get_exp_dat("Cacna1h")
armc4_exp <- get_exp_dat("Armc4")
dscam_exp <- get_exp_dat("Dscam")
spock3_exp <- get_exp_dat("Spock3")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
hhex_exp <- get_exp_dat("Hhex")
sst_exp <- get_exp_dat("Sst")
slc16a7_exp <- get_exp_dat("Slc16a7")

ghrelin_corr <- cbind(ghrelin_corr, arg1_exp, ptprz1_exp, kcnd2_exp, cacna1h_exp, armc4_exp,
                      dscam_exp, spock3_exp, nhs_exp, crhr2_exp, hhex_exp, sst_exp, slc16a7_exp)

pcor <- cor(ghrelin_corr, use= "complete.obs")
round(pcor, digits = 2)
corrplot(pcor, order = "hclust")

# efnb3 and ghsr have very similar correlations patterns, but none of the chromosome 18
# candidate genes are strongly correlated with ghsr or efnb3! 
