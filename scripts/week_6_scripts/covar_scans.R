# Running QTL scans with covariates
# Sadie Allen
# July 11, 2017
# In this script I will attempt to run QTL scans with covariates

# Try to get it to work with 1, then write a loop 

# Clear environment
rm(list = ls())

##### Load in all necessary libraries, data, functions, etc #####
# Load libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)
library(dbplyr)

# Load data
# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load functions
source("scripts/functions.R")

# GET LIST OF MEDIATOR GENES
sst_exp <- get_exp_dat("Sst")
cacna1h_exp <- get_exp_dat("Cacna1h")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
kcnd2_exp <- get_exp_dat("Sst")
arg1_exp <- get_exp_dat("Arg1")
mediator_list <- cbind(ghrelin_list, sst_exp, cacna1h_exp, nhs_exp, crhr2_exp, kcnd2_exp,
                      arg1_exp)

# Rerun important QTL scans
# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ghsr, map = map, threshold = 6, drop = 1.5)
# Peak: chr 18, pos 4.474313, ci: 0.026149-10.21129
# Peak: chr 4, 11.984686, ci: 0.207518-13.32408

qtl.food_ave <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "food_ave", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_ave, map = map, threshold = 6, drop = 1.5)
# Peak: chr 1, pos 172.2283, ci 171.3744-172.4475
# Peak: chr 7, pos 136.4545, ci = 135.5631-136.9891

qtl.efnb3 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.efnb3, map = map, threshold = 6, drop = 1.5)
# chr 4, pos 5.566854, ci 0.0207518 - 13.32408
# chr 11, pos 19.23, ci 17.680768 - 19.61891
# chr 18, pos 3.351, ci 0.026149 - 10.25686


# Create covariates
temp = merge(annot.samples, mediator_list, by = "row.names")
rownames(temp) = temp[,1]
# Now, create the covariates

add_covar_ghsr <- model.matrix(~Sex + Generation + diet_days + ghsr_exp, data = temp)[,-1]
add_covar_sst <- model.matrix(~Sex + Generation + diet_days + sst_exp, data = temp)[,-1]
add_covar_cacna1h <- model.matrix(~Sex + Generation + diet_days + cacna1h_exp, data = temp)[,-1]
add_covar_nhs <- model.matrix(~Sex + Generation + diet_days + nhs_exp, data = temp)[,-1]
add_covar_crhr2 <- model.matrix(~Sex + Generation + diet_days + crhr2_exp, data = temp)[,-1]
add_covar_ghsr <- model.matrix(~Sex + Generation + diet_days + ghsr_exp, data = temp)[,-1]
add_covar_kcnd2 <- model.matrix(~Sex + Generation + diet_days + crhr2_exp, data = temp)[,-1]
add_covar_efnb3 <- model.matrix(~Sex + Generation + diet_days + efnb3_exp, data = temp)[,-1]
add_covar_arg1 <- model.matrix(~Sex + Generation + diet_days + arg1_exp, data = temp)[,-1]

# Mediators of Efnb3: Ghsr, sst, cacna1h, Nhs, Crhr2
qtl.efnb3_ghsrcovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar_ghsr, cores = 4)
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3 vs ghsrcovar", col = "black")
plot_scan1(x = qtl.efnb3_ghsrcovar, map = map, add = TRUE, col = "red")

qtl.efnb3_sstcovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                            kinship = kin, addcovar = add_covar_sst, cores = 4)
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3 vs sstcovar", col = "black")
plot_scan1(x = qtl.efnb3_sstcovar, map = map, add = TRUE, col = "red")

qtl.efnb3_cacna1hcovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                            kinship = kin, addcovar = add_covar_cacna1h, cores = 4)
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3 vs cacna1hcovar", col = "black")
plot_scan1(x = qtl.efnb3_cacna1hcovar, map = map, add = TRUE, col = "red")

qtl.efnb3_nhscovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                                kinship = kin, addcovar = add_covar_nhs, cores = 4)
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3 vs nhscovar", col = "black")
plot_scan1(x = qtl.efnb3_nhscovar, map = map, add = TRUE, col = "red")

qtl.efnb3_crhr2covar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                                kinship = kin, addcovar = add_covar_crhr2, cores = 4)
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3 vs crhr2covar", col = "black")
plot_scan1(x = qtl.efnb3_crhr2covar, map = map, add = TRUE, col = "red")

# Mediators of Ghsr: nhs, cacna1h, kcnd2, efnb3, arg1
qtl.ghsr_nhscovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                           kinship = kin, addcovar = add_covar_nhs, cores = 4)
plot_scan1(x = qtl.ghsr, map = map, main = "ghsr vs nhscovar", col = "black")
plot_scan1(x = qtl.ghsr_nhscovar, map = map, add = TRUE, col = "red")

qtl.ghsr_cacna1hcovar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                           kinship = kin, addcovar = add_covar_cacna1h, cores = 4)
plot_scan1(x = qtl.ghsr, map = map, main = "ghsr vs cacna1hcovar", col = "black")
plot_scan1(x = qtl.ghsr_cacna1hcovar, map = map, add = TRUE, col = "red")

qtl.ghsr_kcnd2covar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                           kinship = kin, addcovar = add_covar_kcnd2, cores = 4)
plot_scan1(x = qtl.ghsr, map = map, main = "ghsr vs kcnd2covar", col = "black")
plot_scan1(x = qtl.ghsr_kcnd2covar, map = map, add = TRUE, col = "red")

qtl.ghsr_efnb3covar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                           kinship = kin, addcovar = add_covar_efnb3, cores = 4)
plot_scan1(x = qtl.ghsr, map = map, main = "ghsr vs efnb3covar", col = "black")
plot_scan1(x = qtl.ghsr_efnb3covar, map = map, add = TRUE, col = "red")

qtl.ghsr_arg1covar <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                           kinship = kin, addcovar = add_covar_arg1, cores = 4)
plot_scan1(x = qtl.ghsr, map = map, main = "ghsr vs arg1covar", col = "black")
plot_scan1(x = qtl.ghsr_arg1covar, map = map, add = TRUE, col = "red")



# Yayyy got it to work for 1! Make a loop
# Goal: loop through genome to see which genes drop a specific lod peak
# Wait... thats gonna take like three million years... how about noooo (the web
# eQTL viewer already does this for you!)
# Ghsr chromosome 18 lod peak: position 4.474313, lod 7.822138
gene_names <- get_gene_names()

for(i in 1:length(gene_names)) {
  gene <- get_exp_dat(gene_names[i])
}

















