# Effect Plots 
# Sadie Allen
# July 5, 2017
# In this script I will create effect plots for 
# phenotypes and gene expressions I have deemed significant to my project
# phenotypes: food_ave, weight_sac, weight_6wk, weight_change
# gene expressions: Ghsr, Svil, Fzd8, Mpp7, Ccny, Mtpap, Usp14, Efnb3

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

# Prepare data
weight_6wk <- matched_phenos$weight_6wk
weight_sac <- matched_phenos$weight_sac
weight_change <- weight_sac - weight_6wk
food_ave <- matched_phenos$food_ave

ghsr_exp <- get_exp_dat("Ghsr")
svil_exp <- get_exp_dat("Svil")
fzd8_exp <- get_exp_dat("Fzd8")
mpp7_exp <- get_exp_dat("Mpp7")
ccny_exp <- get_exp_dat("Ccny")
mtpap_exp <- get_exp_dat("Mtpap")
usp14_exp <- get_exp_dat("Usp14")
efnb3_exp <- get_exp_dat("Efnb3")

phenos <- data.frame(sex = sex, mouse_id = mouse_id, weight_6wk = weight_6wk, weight_sac = weight_sac, weight_change = weight_change,
                     food_ave = food_ave, ghsr = ghsr_exp, svil = svil_exp, fzd8 = fzd8_exp,
                     mpp7 = mpp7_exp, ccny = ccny_exp, mtpap = mtpap_exp, usp14 = usp14_exp, 
                     efnb3 = efnb3_exp)
rownames(phenos) <- phenos$mouse_id

# Prepare to run QTL scans

load("data/qtl_prep.RData")

# Convert genoprobs from DOQTL to QTL2 format
#probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

qtl.weight_6wk <- scan1(genoprobs = probs, pheno = phenos[,3, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_sac <- scan1(genoprobs = probs, pheno = phenos[,4, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.weight_change <- scan1(genoprobs = probs, pheno = phenos[,5, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.food_ave <- scan1(genoprobs = probs, pheno = phenos[,6, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.ghsr <- scan1(genoprobs = probs, pheno = phenos[,7, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.svil <- scan1(genoprobs = probs, pheno = phenos[,8, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.fzd8 <- scan1(genoprobs = probs, pheno = phenos[,9, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.mpp7 <- scan1(genoprobs = probs, pheno = phenos[,10, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.ccny <- scan1(genoprobs = probs, pheno = phenos[,11, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.mtpap <- scan1(genoprobs = probs, pheno = phenos[,12, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.usp14 <- scan1(genoprobs = probs, pheno = phenos[,13, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
qtl.efnb3 <- scan1(genoprobs = probs, pheno = phenos[,14, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)

plot_scan1(x = qtl.weight_6wk, map = map, main = colnames(qtl.weight_6wk)[1])
# chr2
plot_scan1(x = qtl.weight_sac, map = map, main = colnames(qtl.weight_sac)[1])
#chr17
plot_scan1(x = qtl.weight_change, map = map, main = colnames(qtl.weight_change)[1])
#chr8
plot_scan1(x = qtl.food_ave, map = map, main = colnames(qtl.food_ave)[1])
#chr7
plot_scan1(x = qtl.ghsr, map = map, main = colnames(qtl.ghsr)[1])
#chr18
plot_scan1(x = qtl.svil, map = map, main = colnames(qtl.svil)[1])
#chr18
plot_scan1(x = qtl.fzd8, map = map, main = colnames(qtl.fzd8)[1])
# chr18
plot_scan1(x = qtl.mpp7, map = map, main = colnames(qtl.mpp7)[1])
# chr18
plot_scan1(x = qtl.ccny, map = map, main = colnames(qtl.ccny)[1])
# chr18
plot_scan1(x = qtl.mtpap, map = map, main = colnames(qtl.mtpap)[1])
# chr18
plot_scan1(x = qtl.usp14, map = map, main = colnames(qtl.usp14)[1])
# chr18
plot_scan1(x = qtl.efnb3, map = map, main = colnames(qtl.efnb3)[1])
# chr4, chr18

# Effect plots
chr = 2
qtl.weight_6wk.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,3, drop = FALSE],
                                 kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_6wk.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_6wk)[1], scan1_output = qtl.weight_6wk)

qtl.weight_6wk.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos[,3, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_6wk.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_6wk)[1], scan1_output = qtl.weight_6wk)

chr = 17
qtl.weight_sac.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,4, drop = FALSE],
                                 kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_sac.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_sac)[1], scan1_output = qtl.weight_sac)

qtl.weight_sac.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,4, drop = FALSE],
                                kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_sac.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_sac)[1], scan1_output = qtl.weight_sac)

chr = 8
qtl.weight_change.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,5, drop = FALSE],
                                    kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_change.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_change)[1], scan1_output = qtl.weight_change)

qtl.weight_change.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,5, drop = FALSE],
                                kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.weight_change.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weight_change)[1], scan1_output = qtl.weight_change)

chr = 7
qtl.food_ave.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,6, drop = FALSE],
                                    kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.food_ave.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.food_ave)[1], scan1_output = qtl.food_ave)

qtl.food_ave.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,6, drop = FALSE],
                                   kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.food_ave.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.food_ave)[1], scan1_output = qtl.food_ave)

chr = 18
qtl.ghsr.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,7, drop = FALSE],
                               kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.ghsr.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ghsr)[1], scan1_output = qtl.ghsr)

qtl.ghsr.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,7, drop = FALSE],
                              kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.ghsr.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ghsr)[1], scan1_output = qtl.ghsr)

chr = 18
qtl.svil.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,8, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.svil.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.svil)[1], scan1_output = qtl.svil)

qtl.svil.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,8, drop = FALSE],
                          kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.svil.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.svil)[1], scan1_output = qtl.svil)

chr = 18
qtl.fzd8.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,9, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.fzd8.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.fzd8)[1], scan1_output = qtl.fzd8)

qtl.fzd8.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,9, drop = FALSE],
                          kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.fzd8.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.fzd8)[1], scan1_output = qtl.fzd8)

chr = 18
qtl.mpp7.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,10, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.mpp7.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.mpp7)[1], scan1_output = qtl.mpp7)

qtl.mpp7.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,10, drop = FALSE],
                          kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.mpp7.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.mpp7)[1], scan1_output = qtl.mpp7)

chr = 18
qtl.ccny.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,11, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.ccny.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ccny)[1], scan1_output = qtl.ccny)

qtl.ccny.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,11, drop = FALSE],
                          kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.ccny.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ccny)[1], scan1_output = qtl.ccny)

chr = 18
qtl.mtpap.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,12, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.mtpap.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.mtpap)[1], scan1_output = qtl.mtpap)

qtl.mtpap.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,12, drop = FALSE],
                          kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.mtpap.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.mtpap)[1], scan1_output = qtl.mtpap)

chr = 18
qtl.usp14.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,13, drop = FALSE],
                            kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.usp14.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.usp14)[1], scan1_output = qtl.usp14)

qtl.usp14.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,13, drop = FALSE],
                           kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.usp14.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.usp14)[1], scan1_output = qtl.usp14)

chr = 18 # REPEAT W/ CHR 4!!
qtl.efnb3.coef <- scan1coef(genoprobs = probs[,chr], pheno = phenos[,14, drop = FALSE],
                            kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.efnb3.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.efnb3)[1], scan1_output = qtl.efnb3)

qtl.efnb3.blup = scan1blup(genoprobs = probs[,chr], pheno = phenos[,14, drop = FALSE],
                           kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.efnb3.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.efnb3)[1], scan1_output = qtl.efnb3)

# phenotypes: food_ave, weight_sac, weight_6wk, weight_change
# gene expressions: Ghsr, Svil, Fzd8, Mpp7, Ccny, Mtpap, Usp14, Efnb3











































