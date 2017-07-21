# Beta Cell Markers Peak Analysis
# Sadie Allen
# July 21, 2017
# Looking at beta cell marker genes to see if there are any shared peaks

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


# Beta Cells:
Ins1_exp <- get_exp_dat("Ins1")
Ins2_exp <- get_exp_dat("Ins2")
Pdx1_exp <- get_exp_dat("Pdx1")
Mafa_exp <- get_exp_dat("Mafa")
Nkx61_exp <- get_exp_dat("Nkx6-1")
Glp1r_exp <- get_exp_dat("Glp1r")
Ucn3_exp <- get_exp_dat("Ucn3")
Slc2a2_exp <- get_exp_dat("Slc2a2")

beta_genes <- cbind(Ins1_exp, Ins2_exp, Pdx1_exp, Mafa_exp, Nkx61_exp,
                    Glp1r_exp, Ucn3_exp, Slc2a2_exp)

# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# QTLs of beta cell markers
qtl.pdx1 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Pdx1_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.glp1r <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Glp1r_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.ins1 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Ins1_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.ins2 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Ins2_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.mafa <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Mafa_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.nkx61 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Nkx61_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.ucn3 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Ucn3_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.slc2a2 <- scan1(genoprobs = probs, pheno = beta_genes[, colnames(beta_genes) == "Slc2a2_exp", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)

quartz()
par(mfrow = c(4, 2))
plot_scan1(x = qtl.pdx1, map = map, main = "pdx1", col = "black")
plot_scan1(x = qtl.glp1r, map = map, main = "glp1r", col = "black")
plot_scan1(x = qtl.ins1, map = map, main = "ins1", col = "black")
plot_scan1(x = qtl.ins2, map = map, main = "ins2", col = "black")
plot_scan1(x = qtl.mafa, map = map, main = "mafa", col = "black")
plot_scan1(x = qtl.nkx61, map = map, main = "nkx61", col = "black")
plot_scan1(x = qtl.ucn3, map = map, main = "ucn3", col = "black")
plot_scan1(x = qtl.slc2a2, map = map, main = "slc2a2", col = "black")
# No common peaks among the genes





















