# Estimating Heritability
# Sadie Allen
# July 7, 2017
# In this script I will determine the heritability of my phenotypes of interest:
# food_ave, weight_sac, weight_6wk, and weight_change

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

# Load data
# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load functions
source("scripts/functions.R")

# QTL scan prep
load("data/qtl_prep.RData")
#probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# food_ave
est_herit(pheno = ghrelin_list[,5], kinship = kin, addcovar = add_covar, cores = 4)
# WHY DIS NO WORK?























