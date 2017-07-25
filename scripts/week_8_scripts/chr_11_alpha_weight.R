# Chromosome 11 Alpha Eigengene and Body Weight Peak
# Sadie Allen
# July 24, 2017
# Determining genes under chromosome 11 peak 

rm(list = ls())

# Load libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)

# Load all functions
source("scripts/functions.R")

# Load data
# islet data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# UPDATED PHENOTYPE DATA
matched_phenos <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv", as.is=TRUE)
rownames(matched_phenos) <- matched_phenos[,1]

########################################
# prepare phenotypes
food_6wk <- matched_phenos$food_6wk
food_ave <- matched_phenos$food_ave
weight_6wk <- matched_phenos$weight_6wk
weight_16wk <- matched_phenos$weight_16wk
weight_change <- weight_16wk - weight_6wk

phenos <- data.frame(food_6wk = food_6wk, food_ave = food_ave, weight_6wk = weight_6wk, 
                     weight_16wk = weight_16wk, weight_change = weight_change)

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

#food_6wk, food_ave, weight_16wk, and weight_6wk appear signifcantly more normal after
# log transformation

phenos$food_6wk <- log10(phenos$food_6wk)
phenos$food_ave <- log10(phenos$food_ave)
phenos$weight_6wk <- log10(phenos$weight_6wk)
phenos$weight_16wk <- log10(phenos$weight_16wk)
phenos$weight_change <- phenos$weight_16wk - phenos$weight_6wk

food_6wk <- matched_phenos$food_6wk
food_ave <- matched_phenos$food_ave
weight_6wk <- matched_phenos$weight_6wk
weight_16wk <- matched_phenos$weight_16wk
weight_change <- phenos$weight_change

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

# Yay they are all normal!!! :))

########################################

alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

gene_phenos <- cbind(phenos, alpha_eigengene)
rownames(gene_phenos) <- annot.samples$Mouse.ID

####################################
peak_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/results/islet_composition/alpha_11_genes.csv")

peak_genes_pc <- peak_genes[peak_genes$Gene.type == "protein_coding",]
peak_genes_pc

alpha_data <- read.csv(file = "data/islet_composition/alpha_module.csv")
alpha_genes <- alpha_data$gene_symbol

intersect(peak_genes_pc$Gene.name, alpha_genes)
# No shared genes! Good... because I might be getting at the driving gene??
###################################

load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# run initial QTL: 
qtl.alpha <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)

## Get gene expressions
vstm2a_exp <- get_exp_dat("Vstm2a")
sec61g_exp <- get_exp_dat("Sec61g")
egfr_exp <- get_exp_dat("Egfr")
plek_exp <- get_exp_dat("Plek")
cnrip1_exp <- get_exp_dat("Cnrip1")
pno1_exp <- get_exp_dat("Pno1")
c1d_exp <- get_exp_dat("C1d")
etaa1_exp <- get_exp_dat("Etaa1")
meis1_exp <- get_exp_dat("Meis1")

gene_list <- cbind(vstm2a_exp, sec61g_exp, egfr_exp, plek_exp, cnrip1_exp, pno1_exp, c1d_exp,
                   etaa1_exp, meis1_exp)
rownames(gene_list) <- annot.samples$Mouse.ID

# create new covariates:
diet_days <- matched_phenos$diet_days
temp = merge(annot.samples, gene_list, by = "row.names")
rownames(temp) = temp[,1]
# Now, create the covariates

add_covar_vstm2a <- model.matrix(~Sex + Generation + diet_days + vstm2a_exp, data = temp)[,-1]
add_covar_sec61g <- model.matrix(~Sex + Generation + diet_days + sec61g_exp, data = temp)[,-1]
add_covar_egfr <- model.matrix(~Sex + Generation + diet_days + egfr_exp, data = temp)[,-1]
add_covar_plek <- model.matrix(~Sex + Generation + diet_days + plek_exp, data = temp)[,-1]
add_covar_cnrip1 <- model.matrix(~Sex + Generation + diet_days + cnrip1_exp, data = temp)[,-1]
add_covar_pno1 <- model.matrix(~Sex + Generation + diet_days + pno1_exp, data = temp)[,-1]
add_covar_c1d <- model.matrix(~Sex + Generation + diet_days + c1d_exp, data = temp)[,-1]
add_covar_etaa1 <- model.matrix(~Sex + Generation + diet_days + etaa1_exp, data = temp)[,-1]
add_covar_meis1 <- model.matrix(~Sex + Generation + diet_days + meis1_exp, data = temp)[,-1]

# Covariate Scans

qtl.alpha_vstm2acovar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                             kinship = kin, addcovar = add_covar_vstm2a, cores = 4)
qtl.alpha_sec61gcovar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_sec61g, cores = 4)
qtl.alpha_egfrcovar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_egfr, cores = 4)
qtl.alpha_plekcovar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_plek, cores = 4)
qtl.alpha_cnrip1covar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_cnrip1, cores = 4)
qtl.alpha_pno1covar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_pno1, cores = 4)
qtl.alpha_c1dcovar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_c1d, cores = 4)
qtl.alpha_etaa1covar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_etaa1, cores = 4)
qtl.alpha_meis1covar <- scan1(genoprobs = probs, pheno = gene_phenos[, colnames(gene_phenos) == "alpha_eigengene", drop = FALSE],
                               kinship = kin, addcovar = add_covar_meis1, cores = 4)

head(qtl.alpha)
qtl.alpha_11 <- qtl.alpha[rownames(qtl.alpha) == grep("^11_", rownames(qtl.alpha))]

rownames_11 <- grep("^11_", rownames(qtl.alpha))

qtl.alpha_11 <- qtl.alpha[rownames_11, ]
max(qtl.alpha_11)
which(qtl.alpha_11 == max(qtl.alpha_11))

# rowname for max lod: 11_16653937 
# THIS IS WHERE I CHECK PEAK SIZE

qtl.alpha["11_16653937",]
# 6.51551
qtl.alpha_vstm2acovar["11_16653937",]
# 6.515146
qtl.alpha_sec61gcovar["11_16653937",]
# 6.256009
qtl.alpha_egfrcovar["11_16653937",]
# 6.352032
qtl.alpha_plekcovar["11_16653937",]
# 6.531509
qtl.alpha_cnrip1covar["11_16653937",]
# 6.670112
qtl.alpha_pno1covar["11_16653937",]
# 3.47274  PEAK DROP!!
qtl.alpha_c1dcovar["11_16653937",]
# 8.360647
qtl.alpha_etaa1covar["11_16653937",]
# 6.84806
qtl.alpha_meis1covar["11_16653937",]
# 6.492837

# One gene dropped the peak! Pno1
# relationship to confirm: QA11 --> pno1_exp --> alpha_eigengene
sex <- annot.samples$Sex
QA11 <- get_genoprob(11, 16.65394)

# alpha eigengene linked to QA11
anova(lm(alpha_eigengene~sex + QA11)) # 0.0002085 ***
# alpha eigengene is linked to Pno1 expression... i hope
anova(lm(alpha_eigengene~sex + pno1_exp)) # 7.383e-09 *** 
# alpha eigengene not linked to QA11 after accounting for pno1_exp
anova(lm(alpha_eigengene~sex + pno1_exp + QA11)) # 0.0128*
# alpha eigengene still linked to pno1_exp after accounting for QA11
anova(lm(alpha_eigengene~sex  + QA11 + pno1_exp)) # 1.936e-06 ***
# This looks pretty good 

quartz()
plot_scan1(x = qtl.alpha, map = map, main = "alpha eigengene with Pno1 covar", col = "black")
plot_scan1(x = qtl.alpha_pno1covar, map = map, add = TRUE, col = "red")


















