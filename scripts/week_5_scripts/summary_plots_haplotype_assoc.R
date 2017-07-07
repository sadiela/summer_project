# Summary
# Sadie Allen
# July 6, 2017
# QTL scans of genes and phenotypes of interest for meeting with Gary 

# clin pheno scans
# gene expression traits
# effect plots for qtl peaks

# Clear environment
rm(list = ls())

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
#matched_phenos <- read.csv("data/matched_phenos.csv")
#rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load functions
source("scripts/functions.R")

#colnames(ghrelin_list)
#ghrelin_list$X <- NULL
#ghrelin_list$Glu_0min <- NULL
#ghrelin_list$Glu_sac <- NULL
#ghrelin_list$Ins_0min <- NULL
#ghrelin_list$Ins_sac <- NULL

# Get additional phenotypes
#svil_exp <- get_exp_dat("Svil")
#fzd8_exp <- get_exp_dat("Fzd8")
#mpp7_exp <- get_exp_dat("Mpp7")
#ccny_exp <- get_exp_dat("Ccny")
#mtpap_exp <- get_exp_dat("Mtpap")
#efnb3_exp <- get_exp_dat("Efnb3")
#weight_6wk <- ghrelin_list$weight_6wk
#weight_change <- ghrelin_list$weight_sac - ghrelin_list$weight_6wk
#ghrelin_list <- cbind(ghrelin_list, weight_6wk, weight_change, svil_exp, fzd8_exp, mpp7_exp, ccny_exp, mtpap_exp, efnb3_exp)
#names(ghrelin_list)
#write.csv(x = ghrelin_list, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")

# QTL scan prep
load("data/qtl_prep.RData")
#probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

names(ghrelin_list)

qtl.food_ave <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "food_ave", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_ave, map = map, threshold = 6, drop = 1.5)

qtl.weight_sac <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_sac", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_sac, map = map, threshold = 6, drop = 1.5)

qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ghsr, map = map, threshold = 6, drop = 1.5)

qtl.weight_6wk <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_6wk", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_6wk, map = map, threshold = 6, drop = 1.5)

qtl.weight_change <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_change", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_change, map = map, threshold = 6, drop = 1.5)

qtl.svil <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "svil_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.svil, map = map, threshold = 6, drop = 1.5)

qtl.fzd8 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "fzd8_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.fzd8, map = map, threshold = 6, drop = 1.5)

qtl.mpp7 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "mpp7_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.mpp7, map = map, threshold = 6, drop = 1.5)

qtl.ccny <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ccny_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ccny, map = map, threshold = 6, drop = 1.5)

qtl.mtpap <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "mtpap_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.mtpap, map = map, threshold = 6, drop = 1.5)

qtl.efnb3 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.efnb3, map = map, threshold = 6, drop = 1.5)



quartz()
par(mfrow = c(6, 2))

plot_scan1(x = qtl.food_ave, map = map, main=colnames(qtl.food_ave)[1])
plot_scan1(x = qtl.weight_sac, map = map, main=colnames(qtl.weight_sac)[1])
plot_scan1(x = qtl.ghsr, map = map, main=colnames(qtl.ghsr)[1])
plot_scan1(x = qtl.weight_6wk, map = map, main=colnames(qtl.weight_6wk)[1])
plot_scan1(x = qtl.weight_change, map = map, main=colnames(qtl.weight_change)[1])
plot_scan1(x = qtl.svil, map = map, main=colnames(qtl.svil)[1])
plot_scan1(x = qtl.fzd8, map = map, main=colnames(qtl.fzd8)[1])
plot_scan1(x = qtl.mpp7, map = map, main=colnames(qtl.mpp7)[1])
plot_scan1(x = qtl.ccny, map = map, main=colnames(qtl.ccny)[1])
plot_scan1(x = qtl.mtpap, map = map, main=colnames(qtl.mtpap)[1])
plot_scan1(x = qtl.efnb3, map = map, main=colnames(qtl.efnb3)[1])


# Haplotype association for peaks

# Add gene plot
# create data structure for gene plot
# this will allow you to have the genes under the haplotype association
genes <- data.frame(chr = annot.mrna$chr, 
                    start = annot.mrna$start,
                    stop = annot.mrna$end, 
                    strand = annot.mrna$strand,
                    Name = annot.mrna$symbol, 
                    stringsAsFactors = FALSE)

# Weight_sac:
chr = 17
assoc.weight_sac <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 6, 
                                  addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                  start = 31.318, end = 78.823, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.weight_sac[[1]], snpinfo = assoc.weight_sac[[2]], 
             drop.hilit = 1, xlim = c(31.318, 36))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 31.318e6 & genes$stop < 36e6,], 
           xlim = c(31.318, 36))

# Ghsr
chr = 18
assoc.ghsr_18 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 7, 
                                  addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                  start = 0.2075, end = 13.3241, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.ghsr_18[[1]], snpinfo = assoc.ghsr_18[[2]], 
             drop.hilit = 1, xlim = c(0.2075, 13.3241))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 0.2075e6 & genes$stop < 13.3241e6,], 
           xlim = c(0.2075, 13.3241))

chr = 4
assoc.ghsr_4 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 7, 
                                  addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                  start = 0.026, end = 10.22, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.ghsr_4[[1]], snpinfo = assoc.ghsr_4[[2]], 
             drop.hilit = 1, xlim = c(0.026, 10.22))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 0.026e6 & genes$stop < 10.22e6,], 
           xlim = c(0.026, 10.22))

chr = 1
assoc.weight_6wk_1 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 8, 
                              addcovar = add_covar, k = kin, markers = snps, chr = chr,
                              start = 158, end = 160, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.weight_6wk_1[[1]], snpinfo = assoc.weight_6wk_1[[2]], 
             drop.hilit = 1, xlim = c(158, 160))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 158e6 & genes$stop < 160e6,], 
           xlim = c(158, 160))

chr = 11
assoc.weight_6wk_11 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 8, 
                                    addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                    start = 10.28, end = 19.5143, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.weight_6wk_11[[1]], snpinfo = assoc.weight_6wk_11[[2]], 
             drop.hilit = 1, xlim = c(10.28, 19.5143))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 10.28e6 & genes$stop < 19.5143e6,], 
           xlim = c(10.28, 19.5143))

chr = 8
assoc.weight_change_8 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 9, 
                                     addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                     start = 109, end = 111, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.weight_change_8[[1]], snpinfo = assoc.weight_change_8[[2]], 
             drop.hilit = 1, xlim = c(109, 111))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 109e6 & genes$stop < 111e6,], 
           xlim = c(109, 111))


# I am not going to do association mapping for local QTLs... I don't think it makes sense to do so
chr = 18
assoc.efnb3_18 <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 15, 
                                       addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                       start = .026, end = 10.3, ncores = 4)
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.efnb3_18[[1]], snpinfo = assoc.efnb3_18[[2]], 
             drop.hilit = 1, xlim = c(.026, 10.3))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > .026e6 & genes$stop < 10.3e6,], 
           xlim = c(.026, 10.3))



















