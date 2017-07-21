# Gene cluster confirmation
# Sadie Allen
# July 21, 2017
# Checking that expressions of gene markers of alpha, beta, and delta cells form correlation clusters

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

# I will take the top 8 genes from the alpha and delta modules as well as 8 genes known to 
# be beta cell markers and create a correlation matrix to establish grouping
# Delta Cells:
ptprz1_exp <- get_exp_dat("Ptprz1")
ghsr_exp <- get_exp_dat("Ghsr")
hhex_exp <- get_exp_dat("Hhex")
rbp4_exp <- get_exp_dat("Rbp4")
gap43_exp <- get_exp_dat("Gap43")
f5_exp <- get_exp_dat("F5")
nhs_exp <- get_exp_dat("Nhs")
efnb3_exp <- get_exp_dat("Efnb3")

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

alphabetadelta_phenos <- cbind(ptprz1_exp, ghsr_exp, hhex_exp, rbp4_exp, gap43_exp,
                               f5_exp, nhs_exp, efnb3_exp, Irx2_exp, Gcg_exp, Sgce_exp, Sstr2_exp, Mafb_exp,
                               Arx_exp, Galnt13_exp, Mamld1_exp, Ins1_exp, Ins2_exp, Pdx1_exp, Mafa_exp, Nkx61_exp,
                               Glp1r_exp, Ucn3_exp, Slc2a2_exp)

pcor <- cor(alphabetadelta_phenos, use= "complete.obs")
#round(pcor, digits=2)
quartz()
corrplot(pcor, order = "hclust")

############################################
# How do these gene clusters correlate with my phenotypes of interest?
############################################

# Prepare pheno data
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

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

# Yay they are all normal!!! :))

# Load in eigengene data
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

gene_frame <- data.frame(alpha = alpha_eigengene, delta = delta_eigengene, beta = Ins2_exp)

gene_pheno <- cbind(gene_frame, phenos)
rownames(gene_pheno) <- annot.samples$Mouse.ID

pcor <- cor(gene_pheno, use= "complete.obs")
#round(pcor, digits=2)
quartz()
corrplot(pcor, order = "hclust")































