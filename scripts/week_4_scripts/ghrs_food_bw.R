# Exploring the relationship between ghrelin receptors, average food intake, and bodyweight
# Sadie Allen
# Jun 26, 2017
# In this script I will delve deeper into the relationships between ghrelin receptors, 
# average food intake, and bodyweight.

rm(list = ls())

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
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load all functions
source("scripts/functions.R")


## Step 1: Find QTL peak for ghsr (14) ##
# GARY QUESTION: What is the threshold for a significant QTL peak? 6, 8, 11

###### QTL PREP ######## (needed for all qtl scans)
# Convert genoprobs from DOQTL to QTL2 format
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Calculate kinship
kin <- calc_kinship(probs = probs, type = "loco", cores = 4)

# Additive covariates
temp = merge(annot.samples, ghrelin_shortlist, by = "row.names")
rownames(temp) = temp[,1]
# annot.samples <- merge(annot.samples, diet_days_id, by.x="Mouse.ID", by.y = "Mouse.ID")
# Now, create the covariate
add_covar <- model.matrix(~Sex + Generation + diet_days, data = temp)[,-1]

# Convert a marker map organized as data frame to a list
map <- map_df_to_list(map = snps, pos_column = "bp")


### RUNNING AND PLOTTING SCANS ###

# Perform scan
qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_shortlist[,colnames(ghrelin_shortlist)=="ghsr", drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
# Plot scan
plot_scan1(x = qtl.ghsr, map = map, main = colnames(qtl.ghsr)[1])
# Find peak positions
find_peaks(qtl.ghsr, map, threshold = 6, drop = 1.5) # What does the drop argument do?
# chromosome 4, pos 11.984686, lod = 6.859136 (interval 0.207518 - 13.32408)
# chromosome 18, pos 4.474313, lod = 7.822138 (interval 0.26149 - 10.21129)

Nhs_exp <- get_exp_dat("Nhs")
ghrelin_shortlist <- cbind(ghrelin_shortlist, Nhs_exp)
qtl.Nhs <- scan1(genoprobs = probs, pheno = ghrelin_shortlist[,colnames(ghrelin_shortlist)=="Nhs_exp", drop = FALSE],
                 kinship = kin, addcovar = add_covar, cores = 4)
plot_scan1(x = qtl.Nhs, map = map, main = colnames(qtl.Nhs)[1])
find_peaks(qtl.Nhs, map, threshold = 6, drop = 1.5) # What does the drop argument do?
# chromosome 18, pos 43.10387, LOD 9.471 (located on the X chromosome...)


# GARY QUESTION: How do I get the genotypes at the peaks?

## Step 2: Find genes that have a significant effect on Ghsr, food_ave, and weight_sac ## 

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


# Ghsr
# Nhs, Ptprz1, Hhex, Fgf14, Arg1, Efnb3, Slc16a7, Rbp4 (most significant genes)
# food_ave
# St8sia2, Bsdc1, E2f1, Adam1a, Wdr45, Ebna1bp2, Slc46a3, Dnajc11, Prmt1, Tcp11l2
# weight_sac
# Fkbp11, Nucb2, Klhl24, Wdr45, Ormdl3, Ssr4, Ins2, Klf11, Dapl1, Gm16515


# Check number of transcripts for a particular gene (if very low, data isn't trustworthy)
gene_name <- INSERT_DESIRED_NAME
gene_id <- annot.mrna$id[annot.mrna$symbol == gene_name]
expr_data <- raw.mrna[, colnames(raw.mrna) == gene_id]






