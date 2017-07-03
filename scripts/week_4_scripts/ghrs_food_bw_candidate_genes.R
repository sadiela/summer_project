# Exploring the relationship between ghrelin receptors, average food intake, and bodyweight
# Sadie Allen
# Jun 26, 2017
# In this script I will find genes that have a significant effect on Ghsr, food_ave, and weight_sac

# prep to run script by clearing environment and setting the proper working directory
rm(list = ls())
setwd("/Users/s-allens/Documents/ssp/summer_project/")

## Load Libraries ##
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)

## Load Data ##
# Phenotype Data
ghrelin_list <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list[,2]
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load all functions
source("scripts/functions.R")

# Generates a list of the gene names
gene_names <- get_gene_names()

## Ghsr
# Calculate pvalues for all gene expressions
ghsr_pvals <- pvals_changing_covar(ghrelin_shortlist$ghsr, rankz.mrna, ghrelin_shortlist$sex)
names(ghsr_pvals) <- gene_names

# Get rid of insignificant genes and sort by p-value
ghsr_pvals <- sig_list(ghsr_pvals)
head(ghsr_pvals, n = 10)
# Nhs, Ptprz1, Hhex, Fgf14, Arg1, Efnb3, Slc16a7, Rbp4 (most significant genes)

# Get expression data for several significant genes
Nhs_exp <- get_exp_dat("Nhs")
Ptprz1_exp <- get_exp_dat("Ptprz1")
Hhex_exp <- get_exp_dat("Hhex")
Fgf14_exp <- get_exp_dat("Fgf14")
Arg1_exp <- get_exp_dat("Arg1")
Efnb3_exp <- get_exp_dat("Efnb3")
Rbp4_exp <- get_exp_dat("Rbp4")
Slc16a7_exp <- get_exp_dat("Slc16a7")

# add to pheno list
ghrelin_shortlist <- cbind(ghrelin_shortlist, Nhs_exp, Ptprz1_exp, Hhex_exp,
                           Fgf14_exp, Arg1_exp, Efnb3_exp, Rbp4_exp, Slc16a7_exp)

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Nhs_exp, group = sex, fill = sex, col = sex)) +
                        geom_point() + theme_classic() + geom_smooth(method = "lm", formula = y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Nhs_exp)
# r = 0.8044818

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Ptprz1_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Ptprz1_exp)
# r = 0.806445

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Hhex_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Hhex_exp)
# r = 0.7790784

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Fgf14_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Fgf14_exp)
# r = 0.7790784

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Arg1_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Arg1_exp)
# r = 0.7602128

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Efnb3_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Efnb3_exp)
# r = 0.7679838

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Rbp4_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Rbp4_exp)
# r = 0.737329

ggplot(ghrelin_shortlist, aes(x = ghsr, y = Slc16a7_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_classic() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Slc16a7_exp)
# r = 0.7450842

# Found some promising candidate genes related to ghsr! :)

## food_ave 

# Load data table
pheno_gene_exp_signif <- read.csv("/Users/s-allens/Documents/ssp/summer_project/results/prelim_research/pheno_gene_exp_signif.csv")
pheno_gene_exp_signif$genes <- as.character(pheno_gene_exp_signif$genes)
pheno_gene_exp_signif$X <- NULL
gene_names <- pheno_gene_exp_signif[,1]

# Food_ave
food_ave_pvals <- pheno_gene_exp_signif[,6]
names(food_ave_pvals) <- gene_names
# Get rid of insignificant genes and sort by p-value
food_ave_pvals <- sig_list(food_ave_pvals)
head(food_ave_pvals, n = 10)
# St8sia2, Bsdc1, E2f1, Adam1a, Wdr45, Ebna1bp2, Slc46a3, Dnajc11, Prmt1, Tcp11l2
St8sia2_exp <- get_exp_dat("St8sia2")
Bsdc1_exp <- get_exp_dat("Bsdc1")
E2f1_exp <- get_exp_dat("E2f1")
Adam1a_exp <- get_exp_dat("Adam1a")
Wdr45_exp <- get_exp_dat("Wdr45")
Ebna1bp2_exp <- get_exp_dat("Ebna1bp2")
Slc46a3_exp <- get_exp_dat("Slc46a3")
Dnsjc11_exp <- get_exp_dat("Dnajc11")

ghrelin_shortlist <- cbind(ghrelin_shortlist, St8sia2_exp, Bsdc1_exp, E2f1_exp,
                           Adam1a_exp, Wdr45_exp, Ebna1bp2_exp, Slc46a3_exp, 
                           Dnsjc11_exp)

cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$St8sia2_exp)
#0.53
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Bsdc1_exp)
#-0.54
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$E2f1_exp)
#-0.54
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Adam1a_exp)
#0.58
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Wdr45_exp)
#-0.54
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Ebna1bp2_exp)
#0.52
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Slc46a3_exp)
#-0.47
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Dnsjc11_exp)
#0.47


## Weight_sac
weight_sac_pvals <- pheno_gene_exp_signif[,7]
names(weight_sac_pvals) <- gene_names
# Get rid of insignificant genes and sort by p-value
weight_sac_pvals <- sig_list(weight_sac_pvals)
head(weight_sac_pvals, n = 10)
# Fkbp11, Nucb2, Klhl24, Wdr45, Ormdl3, Ssr4, Ins2, Klf11, Dapl1, Gm16515
#Wdr45 is shared between food_ave and weight_sac!
Fkbp11_exp <- get_exp_dat("Fkbp11")
Nucb2_exp <- get_exp_dat("Nucb2")
Klhl24_exp <- get_exp_dat("Klhl24")
Wdr45_exp <- get_exp_dat("Wdr45")
Ormdl3_exp <- get_exp_dat("Ormdl3")
Ssr4_exp <- get_exp_dat("Ssr4")
Ins2_exp <- get_exp_dat("Ins2")
Dapl1_exp <- get_exp_dat("Dapl1")

ghrelin_shortlist <- cbind(ghrelin_shortlist, Fkbp11_exp, Nucb2_exp, Klhl24_exp,
                           Ormdl3_exp, Ssr4_exp, Ins2_exp, Dapl1_exp)

cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Fkbp11_exp)
# -0.50
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Nucb2_exp)
# -0.4387
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Klhl24_exp)
# -0.455
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Ormdl3_exp)
# -0.495
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Ssr4_exp)
# -0.438
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Ins2_exp)
# -0.545
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$Dapl1_exp)
# -0.479

# Make sure that all genes are expressed at a significant level to be used accurately in research
# Ghsr
# Nhs, Ptprz1, Hhex, Fgf14, Arg1, Efnb3, Slc16a7, Rbp4 (most significant genes)
num_expressions("Nhs")
num_expressions("Ptprz1")
num_expressions("Hhex")
num_expressions("Fgf14") # min 3... average 89
num_expressions("Arg1")
num_expressions("Efnb3")
num_expressions("Slc16a7")
num_expressions("Rbp4") # min = 552 this is a good one... I think
# food_ave
# St8sia2, Bsdc1, E2f1, Adam1a, Wdr45, Ebna1bp2, Slc46a3, Dnajc11, Prmt1, Tcp11l2
num_expressions("St8sia2") #min 3
num_expressions("Bsdc1") # min 1350, good
num_expressions("E2f1")
num_expressions("Adam1a")
num_expressions("Wdr45")
num_expressions("Ebna1bp2")
num_expressions("Slc46a3")
num_expressions("Dnajc11")
num_expressions("Prmt1")
num_expressions("Tcp11l2")
# weight_sac
# Fkbp11, Nucb2, Klhl24, Wdr45, Ormdl3, Ssr4, Ins2, Klf11, Dapl1, Gm16515
num_expressions("Fkbp11")
num_expressions("Nucb2")
num_expressions("Klhl24")
num_expressions("Ormdl3")
num_expressions("Ssr4")
num_expressions("Ins2")
num_expressions("Klf11")
num_expressions("Dapl1")
num_expressions("Gm16515")
#All of these are fine

## CANDIDATE GENES ##
# Ghsr_exp: Nhs, Ptprz1, Hhex, Fgf14, Arg1, Efnb3, Slc16a7, Rbp4 
# Food_ave: St8sia2, Bsdc1, E2f1, Adam1a, Wdr45, Ebna1bp2, Slc46a3, Dnajc11, Prmt1, Tcp11l2
# Weight_sac: Fkbp11, Nucb2, Klhl24, Wdr45, Ormdl3, Ssr4, Ins2, Klf11, Dapl1, Gm16515
######################

























