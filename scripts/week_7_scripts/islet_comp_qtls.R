# Islet Composition QTLs
# Sadie Allen
# July 21, 2017
# QTLs involved in new project direction

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

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

# Yay they are all normal!!! :))

########################################

# Load in eigengene data
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

gene_frame <- data.frame(alpha = alpha_eigengene, delta = delta_eigengene, beta = Ins2_exp)

gene_pheno <- cbind(gene_frame, phenos)
rownames(gene_pheno) <- annot.samples$Mouse.ID

# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

qtl.alpha <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "alpha", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.alpha, map = map, threshold = 6, drop = 1.5)

qtl.delta <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "delta", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.delta, map = map, threshold = 6, drop = 1.5)

qtl.beta <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "beta", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.beta, map = map, threshold = 6, drop = 1.5)

qtl.food_6wk <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "food_6wk", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_6wk, map = map, threshold = 6, drop = 1.5)

qtl.food_ave <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "food_ave", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_ave, map = map, threshold = 6, drop = 1.5)

qtl.weight_6wk <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "weight_6wk", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_6wk, map = map, threshold = 6, drop = 1.5)

qtl.weight_16wk <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "weight_16wk", drop = FALSE],
                         kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_16wk, map = map, threshold = 6, drop = 1.5)

qtl.weight_change <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "weight_change", drop = FALSE],
                           kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_change, map = map, threshold = 6, drop = 1.5)

quartz()
par(mfrow=c(3,1))
plot_scan1(x = qtl.alpha, map = map, main = "Alpha Eigengene", col = "black")
plot_scan1(x = qtl.beta, map = map, main = "Beta Surrogate (Ins2)", col = "black")
plot_scan1(x = qtl.delta, map = map, main = "Delta Eigengene", col = "black")

# QTLs Driving alpha, beta, and delta gene expression
# Alpha: chr 1 126-136, chr 6 23-27, chr 11 14-20, chr 15 58-67
# Beta: chr 13 70-72, chr20 74-167
# Delta: chr 4 0.2 - 15, chr 6 0.7-16, chr 13 83-98, chr 18 .02-9

# Correlation plots: alpha, beta, and delta gene expressions and clinical phenotypes
sex <- annot.samples$Sex

#####################################
# COVARIATE SCANS


temp = merge(annot.samples, matched_phenos, by = "row.names")
rownames(temp) = temp[,1]
# Now, create the covariates

add_covar_weightsac <- model.matrix(~Sex + Generation + diet_days + weight_sac, data = temp)[,-1]

qtl.alpha_weight <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "alpha", drop = FALSE],
                          kinship = kin, addcovar = add_covar_weightsac, cores = 4)
find_peaks(qtl.alpha_weight, map = map, threshold = 6, drop = 1.5)

qtl.delta_weight <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "delta", drop = FALSE],
                          kinship = kin, addcovar = add_covar_weightsac, cores = 4)

qtl.beta_weight <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "beta", drop = FALSE],
                         kinship = kin, addcovar = add_covar_weightsac, cores = 4)

quartz()
plot_scan1(x = qtl.alpha_weight, map = map, main = "Alpha Eigengene with Weight as Covariate", col = "red")
plot_scan1(x = qtl.alpha, map = map, add = TRUE, col = "black")

quartz()
plot_scan1(x = qtl.delta_weight, map = map, main = "Delta Eigengene with Weight as Covariate", col = "red")
plot_scan1(x = qtl.delta, map = map, add = TRUE, col = "black")

quartz()
plot_scan1(x = qtl.beta, map = map, main = "Ins2 with Weight as Covariate", col = "black")
plot_scan1(x = qtl.beta_weight, map = map, add = TRUE, col = "red")





















