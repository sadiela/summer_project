# Chromosome 11 effect plots: weight phenotypes and alpha eigengene 
# Sadie Allen
# July 24, 2017
# Determining if weight phenos and alpha eigengene have similar effect plots

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

gene_frame <- data.frame(alpha = alpha_eigengene, delta = delta_eigengene)

gene_pheno <- cbind(gene_frame, phenos)
rownames(gene_pheno) <- annot.samples$Mouse.ID

# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"


qtl.alpha <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "alpha", drop = FALSE],
                                kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_6wk <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "weight_6wk", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_16wk <- scan1(genoprobs = probs, pheno = gene_pheno[, colnames(gene_pheno) == "weight_16wk", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)

chr = 11
qtl.alpha.blup <- scan1blup(genoprobs = probs[,chr], pheno = gene_pheno[,1, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
qtl.weight_6wk.blup <- scan1blup(genoprobs = probs[,chr], pheno = gene_pheno[,5, drop = FALSE],
                            kinship = kin[[chr]], addcovar = add_covar)
qtl.weight_16wk.blup <- scan1blup(genoprobs = probs[,chr], pheno = gene_pheno[,6, drop = FALSE],
                            kinship = kin[[chr]], addcovar = add_covar)

quartz()
plot(x = qtl.alpha.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Alpha Eigengene", scan1_output = qtl.alpha)
plot(x = qtl.weight_6wk.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Weight_6wk", scan1_output = qtl.weight_6wk)
plot(x = qtl.weight_16wk.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Weight_16wk", scan1_output = qtl.weight_16wk)


##### Association Mapping ####

assoc.alpha <- assoc_mapping(probs = probs, pheno = gene_pheno, idx = 1, 
                                 addcovar = add_covar, k = kin, markers = snps, chr = 11,
                                 start = 14, end = 20, ncores = 4)
assoc.weight_6wk <- assoc_mapping(probs = probs, pheno = gene_pheno, idx = 5, 
                             addcovar = add_covar, k = kin, markers = snps, chr = 11,
                             start = 14, end = 20, ncores = 4)
assoc.weight_16wk <- assoc_mapping(probs = probs, pheno = gene_pheno, idx = 6, 
                             addcovar = add_covar, k = kin, markers = snps, chr = 11,
                             start = 14, end = 20, ncores = 4)

genes <- data.frame(chr = annot.mrna$chr, 
                    start = annot.mrna$start,
                    stop = annot.mrna$end, 
                    strand = annot.mrna$strand,
                    Name = annot.mrna$symbol, 
                    stringsAsFactors = FALSE)

quartz()
par(mfrow = c(4,1))
plot_snpasso(scan1output = assoc.alpha[[1]], snpinfo = assoc.alpha[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.weight_6wk[[1]], snpinfo = assoc.weight_6wk[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.weight_16wk[[1]], snpinfo = assoc.weight_16wk[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 11 & genes$start > 14e6 & genes$stop < 20e6,], 
           xlim = c(14, 20))

quartz()
plot_scan1(qtl.alpha, map = map, col = "black")













