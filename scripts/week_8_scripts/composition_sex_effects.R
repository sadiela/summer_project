# Sex effects: islet composition
# Sadie Allen
# July 28, 2017

rm(list = ls())

library(tidyverse)
library(plotly)
library(qtl2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)

source("scripts/functions.R")

#Load data
#pheno data
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

##########
# Correlation matrices by sex

# I will take the top 8 genes from the alpha and delta modules as well as 8 genes known to 
# be beta cell markers and create a correlation matrix to establish grouping
# Delta Cells:
Ptprz1_exp <- get_exp_dat("Ptprz1")
Ghsr_exp <- get_exp_dat("Ghsr")
Hhex_exp <- get_exp_dat("Hhex")
Rbp4_exp <- get_exp_dat("Rbp4")
Gap43_exp <- get_exp_dat("Gap43")
F5_exp <- get_exp_dat("F5")
Nhs_exp <- get_exp_dat("Nhs")
Efnb3_exp <- get_exp_dat("Efnb3")

# Alpha Cells:
Irx2_exp <- get_exp_dat("Irx2")
Gcg_exp <- get_exp_dat("Gcg")
Sgce_exp <- get_exp_dat("Sgce")
Sstr2_exp <- get_exp_dat("Sstr2")
Mafb_exp <- get_exp_dat("Mafb")
Arx_exp <- get_exp_dat("Arx")
Galnt13_exp <- get_exp_dat("Galnt13")
Mamld1_exp <- get_exp_dat("Mamld1")

# Beta Cells:
Ins1_exp <- get_exp_dat("Ins1")
Ins2_exp <- get_exp_dat("Ins2")
Pdx1_exp <- get_exp_dat("Pdx1")
Mafa_exp <- get_exp_dat("Mafa")
Nkx61_exp <- get_exp_dat("Nkx6-1")
Glp1r_exp <- get_exp_dat("Glp1r")
Ucn3_exp <- get_exp_dat("Ucn3")
Slc2a2_exp <- get_exp_dat("Slc2a2")


sex <- annot.samples$Sex

alphabetadelta_exp <- data.frame(sex, Pdx1_exp, Glp1r_exp, Ins1_exp, Ins2_exp,  Mafa_exp, Nkx61_exp,
                                  Ucn3_exp, Slc2a2_exp, Mafb_exp, Arx_exp, Irx2_exp, Sstr2_exp, Rbp4_exp, 
                                 Mamld1_exp, Galnt13_exp, Gcg_exp, Sgce_exp, Ptprz1_exp, Hhex_exp,  Ghsr_exp,
                                 Nhs_exp, Efnb3_exp, Gap43_exp, F5_exp)

# Subset data by sex
gender_sep <- split(alphabetadelta_exp, alphabetadelta_exp$sex)
females <- gender_sep$F
males <- gender_sep$M

pcor_f <- cor(females[2:25], use = "complete.obs")
quartz()
corrplot(pcor_f)

pcor_m <- cor(males[2:25], use = "complete.obs")
quartz()
corrplot(pcor_m)

# Load in eigengenes
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

phenos_df <- data.frame(sex, weight_16wk = matched_phenos$weight_16wk, food_ave = matched_phenos$food_ave,
                        alpha_eigengene, delta_eigengene, beta = Ins2_exp)

phenos_gender <- split(phenos_df, phenos_df$sex)
f_phenos <- phenos_gender$F
m_phenos <- phenos_gender$M

pcor_pheno_f <- cor(f_phenos[2:6], use = "complete.obs")
quartz()
corrplot(pcor_pheno_f, order = "hclust")

pcor_pheno_m <- cor(m_phenos[2:6], use = "complete.obs")
quartz()
corrplot(pcor_pheno_m, order = "hclust")


#######################################
# Linear Regression by Sex

# Ins2 expression
anova(lm(alpha_eigengene ~ beta, data = f_phenos))
anova(lm(delta_eigengene ~ beta, data = f_phenos))

anova(lm(alpha_eigengene ~ beta, data = m_phenos))
anova(lm(delta_eigengene ~ beta, data = m_phenos))

# Weight
anova(lm(alpha_eigengene ~ weight_16wk, data = f_phenos))
anova(lm(delta_eigengene ~ weight_16wk, data = f_phenos))

anova(lm(alpha_eigengene ~ weight_16wk, data = m_phenos))
anova(lm(delta_eigengene ~ weight_16wk, data = m_phenos))




#################################
# Run sex separated QTL Scans!
# QTL Prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Run base scans (no extra covariates)
qtl.alpha_female <- scan1(genoprobs = probs, pheno = f_phenos[, colnames(f_phenos) == "alpha_eigengene", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)

qtl.delta_female <- scan1(genoprobs = probs, pheno = f_phenos[, colnames(f_phenos) == "delta_eigengene", drop = FALSE],
                          kinship = kin, addcovar = add_covar, cores = 4)

qtl.alpha_male <- scan1(genoprobs = probs, pheno = m_phenos[, colnames(m_phenos) == "alpha_eigengene", drop = FALSE],
                          kinship = kin, addcovar = add_covar, cores = 4)

qtl.delta_male <- scan1(genoprobs = probs, pheno = m_phenos[, colnames(m_phenos) == "delta_eigengene", drop = FALSE],
                          kinship = kin, addcovar = add_covar, cores = 4)

quartz()
plot_scan1(qtl.alpha_female, map = map, main = "alpha_female", col = "red")
plot_scan1(qtl.alpha_male, map = map, main = "alpha_male", add = TRUE, col = "blue")

plot_scan1(qtl.delta_female, map = map, main = "delta_female", col = "red")
plot_scan1(qtl.delta_male, map = map, main = "delta_male", add = TRUE, col = "blue")
###############################################


QD18 <- get_genoprob(18, 5.990333)
QD6 <- get_genoprob(6, 4.776124)
sex <- annot.samples$Sex
gen <- annot.samples$Generation

anova(lm(delta_eigengene~ sex + gen + QD18 + sex:QD18))
# Q18 is not sex specific

anova(lm(delta_eigengene~ sex + gen + QD6 + sex:QD6))
# Q6 is not sex specific


#################################################
# are the alpha eigengene peaks sex specific?
QA 





















