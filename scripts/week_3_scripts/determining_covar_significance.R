#6/19/2017 (Week3)
#start in project directory
#setwd("/Users/s-allens/Documents/ssp/summer_project")

library(tidyverse)
library(plotly)
library(qtl2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)

#loading in mRNA data 
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")

# getting access functions that are stored in a separate file
source("scripts/functions.R")

#get phenotype data
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)

#which phenotypes/gene expression traits differ between sexes? (p vals for phenotypes)
#graphical analysis
#create a vector of p-values
fit1 <- anova(lm(Glu_10wk ~ sex, data = matched_phenotypes))
# to get p-val, fit1[1,5]
#tells you if the covariate you give it has a significant
#effect on the phenotype

#Goal 1: Get a vector of p values for the effects of gender for 
#each of the ghrelin phenotypes

#phenotypes that will be relevant to my ghrelin study
ghrelin_phenos <- matched_phenotypes %>% select(sex, Glu_0min, Glu_tAUC, Glu_iAUC, Glu_6wk, 
                  Glu_10wk, Glu_14wk, Glu_sac, Ins_0min, Ins_tAUC, Ins_iAUC, Ins_6wk, 
                  Ins_10wk, Ins_14wk, Ins_sac, food_ave, weight_sac, G33_ins_secrete, 
                  G83_ins_secrete, G167_ins_secrete, KCl_G33_ins_secrete, 
                  GLP1_G83_ins_secrete, AA_G83_ins_secrete, PA_G167_ins_secrete, sex)

#pval_vec gets the list of values, significance sorts them by significance
ghrelin_pvals <- pval_vec(ghrelin_phenos, "sex")
sex_insig_ghrelin <- significance(ghrelin_pvals)
#sex has a significant effect on all of these phenotypes (p < 0.05)

#Goal 2: Use P-val functions to test all gene expression data
sex <- matched_phenotypes$sex
mrna_z_data <- data.frame(rankz.mrna)
rankz_and_gender <- cbind(sex, mrna_z_data)

# change all variables to numerical
for(i in 2:ncol(rankz_and_gender)) {
  rankz_and_gender[,i] <- as.numeric(as.character(rankz_and_gender[,i]))
}

#generate one list of genes that are affected and one of those that are not
gene_exp_pvals <- pval_vec(rankz_and_gender, "sex")
sex_signif_genexp <- significance(gene_exp_pvals)
genes_affected <- sex_signif_genexp$significant_effect
genes_unaffected <- sex_signif_genexp$insignificant_effect

#It worked!! :)

#get the actual names of the genes in each list
signif_genes <- character(0)
insig_genes <- character(0)
for(i in 1:length(genes_affected)) {
  gene_name <- annot.mrna$symbol[annot.mrna$id == names(genes_affected[i])]
  signif_genes <- c(signif_genes, gene_name)
  names(genes_affected[i]) <- gene_name
}
for(i in 1:length(genes_unaffected)) {
  gene_name <- annot.mrna$symbol[annot.mrna$id == names(genes_unaffected[i])]
  insig_genes <- c(insig_genes, gene_name)
  names(genes_unaffected[i]) <- gene_name
}
name_p_sig <- data.frame(name = signif_genes, pval = genes_affected)
name_p_insig <- data.frame(name = insig_genes, pval = genes_unaffected)


#Save this data!
#write.csv(name_p_sig, file = "sex_sig_genes")
#write.csv(name_p_insig, file = "sex_insig_genes")

#if you ever need it again:
name_p_sig <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/sex_sig_genes", as.is=TRUE)
name_p_insig <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/sex_insig_genes", as.is=TRUE)

#time to make some histograms?
name_p_insig
name_p_sig
sex_all_expr <- rbind(name_p_insig, name_p_sig)

#long form of dataframe...?
pvals_long <- melt(sex_all_expr$pval, variable.name = "test", value_name = "pval")

#make histogram
hist(pvals_long$value, breaks = 100)
