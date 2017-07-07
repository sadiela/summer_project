# Haplotype Association Mapping
# Sadie Allen
# July 6, 2017
# In this script I will attempt to get association mapping to work for a single phenotype

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
matched_phenos <- read.csv("data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
sex <- matched_phenos$sex
mouse_id <- matched_phenos$Mouse.ID

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

# Run qtl scan
qtl.food_ave <- scan1(genoprobs = probs, pheno = matched_phenos[, colnames(matched_phenos) == "food_ave", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)

# Plot scan
plot_scan1(x = qtl.food_ave, map = map, main=colnames(qtl.food_ave)[1])

find_peaks(qtl.food_ave, map, threshold = 6, drop =1.5)
# chr 7, pos 136.4545, ci: 135.5631, 136.9492

# Effect plot
chr = 7
qtl.food_ave.coef <- scan1coef(genoprobs = probs[,chr], pheno = matched_phenos[, colnames(matched_phenos) == "food_ave", drop = FALSE],
                               kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.food_ave.coef, map = map[[chr]], columns = 1:8, col = CCcolors, 
     main = colnames(qtl.food_ave)[1], scan1_output = qtl.food_ave)

# Blup
qtl.food_ave.blup <- scan1blup(genoprobs = probs[,chr], pheno = matched_phenos[, colnames(matched_phenos) == "food_ave", drop = FALSE],
                              kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.food_ave.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.food_ave)[1], scan1_output = qtl.food_ave)

# Attempting association mapping
# Arguments:
# probs: genoprobs object in qtl2 format (probs)
# pheno: matrix of phenotypes, samples in rows, phenotypes in columns (matched_phenos)
# idx: column index of the phenotype that you want to map (column of food_ave: 55)
# addcovar: covariates matrix as used in scan1(). (add_covar)
# intcovar: covariate to interact with QTL. (N/A)
# k: list of kinship matrices. (kin)
# markers: data.frame containing 4 columns: marker, chr, bp, cM. (snps)
# chr: Chromosome to map on. (chr 7 for food_ave)
# start: start position for mapping in Mb. (????) 135.5631
# end: end position for mapping in Mb. (?????) 136.9492
# ncores: number of cores to use in mapping. (4)
# db.file: Location of the mySQL database containing the Sanger SNPs. ("/Users/s-allens/Documents/ssp/summer_project/data/ccfoundersnps.sqlite")

# FUNCTION:
assoc_mapping = function(probs, pheno, idx, addcovar, intcovar = NULL, k, 
                         markers, chr, start, end, ncores, db.file = "/Users/s-allens/Documents/ssp/summer_project/data/ccfoundersnps.sqlite") {
  # subsetting probs and k
  probs = probs[,chr]
  k = k[[chr]]
  
  # split markers into a vector of map positions
  snps$pos = snps$pos*1e-6
  map = map_df_to_list(map = snps, pos_column = "pos")
  
  # Extract SNPs from the database
  my_db = src_sqlite(db.file, create = FALSE)
  snpinfo = tbl(my_db, sql(paste0("SELECT * FROM snps WHERE chr='", 
                chr, "' AND pos_Mbp>=", start, " AND pos_Mbp<=", end))) %>% collect(n = Inf)
  
  # Replace names for future methods
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")
  
  # Index groups of similar SNPs
  snpinfo = index_snps(map = map, snpinfo) 
  
  # Find which phenotype data exist
  sel = !is.na(pheno[,idx])
  
  # Convert genoprobs to snpprobs
  snppr = genoprob_to_snpprob(probs[sel,], snpinfo)
  
  # Scan1
  assoc = scan1(pheno = pheno[sel, idx, drop = FALSE], kinship = k[sel,sel],
                genoprobs = snppr, addcovar = addcovar[sel,], intcovar = addcovar[sel, intcovar],
                cores = ncores)
  
  # Return scan data
  return(list(assoc, snpinfo))  
}

library(dbplyr)

assoc.food_ave <- assoc_mapping(probs = probs, pheno = matched_phenos, idx = 55, 
                                   addcovar = add_covar, k = kin, markers = snps, chr = 7,
                                   start = 135.5631, end = 136.9492, ncores = 4)

# Association mapping plot.
plot_snpasso(scan1output = assoc.food_ave[[1]], snpinfo = assoc.food_ave[[2]], 
             drop.hilit = 1)

# Association mapping plot with narrower boundaries.
plot_snpasso(scan1output = assoc.food_ave[[1]], snpinfo = assoc.food_ave[[2]], 
             drop.hilit = 1, xlim = c(136.2, 136.7))


# Add gene plot
# create data structure for gene plot
genes <- data.frame(chr = annot.mrna$chr, 
                    start = annot.mrna$start,
                    stop = annot.mrna$end, 
                    strand = annot.mrna$strand,
                    Name = annot.mrna$symbol, 
                    stringsAsFactors = FALSE)

layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))

chr = 7
plot_snpasso(scan1output = assoc.food_ave[[1]], snpinfo = assoc.food_ave[[2]], 
             drop.hilit = 1, xlim = c(134, 138))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 134e6 & genes$stop < 138e6,], 
           xlim = c(134, 138))




# MESS W/ COLORS
#genes.ss = genes[genes$chr == chr & genes$start > 107e6 & genes$stop < 109e6,]
# NOTE: You have to sort the genes by position for this to work.
#genes.ss = genes.ss[order(genes.ss$start),]
#gene.colors = rep("black", nrow(genes.ss))
#gene.colors[genes.ss$Name == "Gnai3"] = "red"

#plot_genes(genes.ss, xlim = c(107, 109), colors = gene.colors)

#plot_genes(genes[genes$chr == chr,], xlim = c(107, 109))



# NOW DO W/ MORE PHENOS
qtl.weight_sac <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_sac", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_sac, threshold = 6, drop = 1.5, map = map)

qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_6wk <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_6wk", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_change <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "weight_change", drop = FALSE],
                           kinship = kin, addcovar = add_covar, cores = 4)
qtl.svil <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "svil_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.fzd8 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "fzd8_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.mpp7 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "mpp7_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.ccny <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ccny_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.mtpap <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "mtpap_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.efnb3 <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "efnb3_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)


assoc.weight_sac <- assoc_mapping(probs = probs, pheno = ghrelin_list, idx = 4, 
                                addcovar = add_covar, k = kin, markers = snps, chr = 17,
                                start = 31.31825, end = 78.82301, ncores = 4)

layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))

chr = 17
plot_snpasso(scan1output = assoc.weight_sac[[1]], snpinfo = assoc.weight_sac[[2]], 
             drop.hilit = 1, xlim = c(31, 79))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 31e6 & genes$stop < 79e6,], 
           xlim = c(31, 79))


