# Back to Delta Cells
# Sadie Allen
# July 24, 2017
# Revisiting peaks on chromosome 6 and 18 for the delta cell eigengene

rm(list = ls())

library(tidyverse)
library(plotly)
library(qtl2)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(corrplot)
library(reshape2)
library(RSQLite)

source("scripts/functions.R")

#Load data
#pheno data
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Run delta cell QTL scans

########################################
# prepare phenotypes
food_6wk <- matched_phenos$food_6wk
food_ave <- matched_phenos$food_ave
weight_6wk <- matched_phenos$weight_6wk
weight_16wk <- matched_phenos$weight_16wk
weight_change <- weight_16wk - weight_6wk

phenos <- data.frame(food_6wk = food_6wk, food_ave = food_ave, weight_6wk = weight_6wk, 
                     weight_16wk = weight_16wk, weight_change = weight_change)

phenos$food_6wk <- log10(phenos$food_6wk)
phenos$food_ave <- log10(phenos$food_ave)
phenos$weight_6wk <- log10(phenos$weight_6wk)
phenos$weight_16wk <- log10(phenos$weight_16wk)
phenos$weight_change <- phenos$weight_16wk - phenos$weight_6wk
########################################
# get data for delta eigengene
delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

# load gene expression data
ptprz1_exp <- get_exp_dat("Ptprz1")
armc4_exp <- get_exp_dat("Armc4")
zfp438_exp <- get_exp_dat("Zfp438")

# add gene expressions to data frame 
phenos_genes <- cbind(phenos, delta_eigengene, ptprz1_exp, armc4_exp, zfp438_exp) 

# QTL Prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Run base scans (no extra covariates)
qtl.delta <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.ptprz1 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "ptprz1_exp", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.armc4 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "armc4_exp", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)
qtl.zfp438 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "zfp438_exp", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)

# find peak locations for delta eigengene 
find_peaks(qtl.delta, map = map, threshold = 6, drop = 1.5)
#   lodindex       lodcolumn chr       pos      lod     ci_lo     ci_hi
# 1        1 delta_eigengene   4 13.186538 6.837581  0.207518 14.424397
# 2        1 delta_eigengene   6  4.776124 9.075374  0.738935 15.498195
# 3        1 delta_eigengene  13 95.992420 6.026294 83.460496 98.259693
# 4        1 delta_eigengene  18  5.990333 8.463098  0.026149  8.530216
#0-16 Mbp

# Plot scans
quartz()
plot_scan1(qtl.delta, map = map, main = "Delta Eigengene", col = "black")
plot_scan1(qtl.ptprz1, map = map, main = "Ptprz1", col = "black")
plot_scan1(qtl.armc4, map = map, main = "Armc4", col = "black")
plot_scan1(qtl.zfp438, map = map, main = "Zfp438", col = "black")

####################################################
# Covariate Scans #
temp = merge(annot.samples, matched_phenos, by = "row.names")
rownames(temp) = temp[,1]

# Now, create the covariates
add_covar_ptprz1 <- model.matrix(~Sex + Generation + diet_days + ptprz1_exp, data = temp)[,-1]
add_covar_armc4 <- model.matrix(~Sex + Generation + diet_days + armc4_exp, data = temp)[,-1]
add_covar_zfp438 <- model.matrix(~Sex + Generation + diet_days + zfp438_exp, data = temp)[,-1]
add_covar_zfp438_armc4 <- model.matrix(~Sex + Generation + diet_days + zfp438_exp + armc4_exp, data = temp)[,-1]

# run scans
qtl.delta_ptprz1 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                   kinship = kin, addcovar = add_covar_ptprz1, cores = 4)
find_peaks(qtl.delta_ptprz1, map = map, threshold = 6, drop = 1.5)
#  lodindex       lodcolumn chr      pos       lod     ci_lo    ci_hi
#1        1 delta_eigengene   6 23.01854 29.571516 22.290482 23.34616 what??!?!
#2        1 delta_eigengene  18 23.82573  6.613041  0.026149 28.13084

qtl.delta_armc4 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                   kinship = kin, addcovar = add_covar_armc4, cores = 4)

qtl.delta_zfp438 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                    kinship = kin, addcovar = add_covar_zfp438, cores = 4)

qtl.delta_zfp438_armc4 <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                          kinship = kin, addcovar = add_covar_zfp438_armc4, cores = 4)

# Plot scans against original
quartz()
plot_scan1(qtl.delta_ptprz1, map = map, main = "Delta Eigengene with Ptprz1 Covariate", col = "red")
plot_scan1(qtl.delta, map = map, add = TRUE, col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.delta, map = map, main = "Delta Eigengene with Armc4 Covariate", col = "black")
plot_scan1(qtl.delta_armc4, map = map, add = TRUE, col = "red")

plot_scan1(qtl.delta, map = map, main = "Delta Eigengene with Zfp438 Covariate", col = "black")
plot_scan1(qtl.delta_zfp438, map = map, add = TRUE, col = "red")

plot_scan1(qtl.delta, map = map, main = "Delta Eigengene with Armc4 and Zfp438 Covariates", col = "black")
plot_scan1(qtl.delta_zfp438_armc4, map = map, add = TRUE, col = "red")

###################
# Effect Plots: Chromosome 6 
chr = 6
qtl.delta_6.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
qtl.delta_ptprz1_6.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                            kinship = kin[[chr]], addcovar = add_covar_ptprz1)

# Effect Plots: Chromosome 18 
chr = 18
qtl.delta_18.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                                     kinship = kin[[chr]], addcovar = add_covar)
qtl.delta_armc4_18.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                                  kinship = kin[[chr]], addcovar = add_covar_armc4)
qtl.delta_zfp438_18.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                                  kinship = kin[[chr]], addcovar = add_covar_zfp438)
qtl.delta_zfp438_armc4_18.blup <- scan1blup(genoprobs = probs[,chr], pheno = phenos_genes[,6, drop = FALSE],
                                  kinship = kin[[chr]], addcovar = add_covar_zfp438_armc4)

# Effect plots: Candidate genes 
qtl.ptprz1.blup <- scan1blup(genoprobs = probs[,6], pheno = phenos_genes[,7, drop = FALSE],
                             kinship = kin[[6]], addcovar = add_covar)
qtl.armc4.blup <- scan1blup(genoprobs = probs[,18], pheno = phenos_genes[,8, drop = FALSE],
                            kinship = kin[[18]], addcovar = add_covar)
qtl.zfp438.blup <- scan1blup(genoprobs = probs[,18], pheno = phenos_genes[,9, drop = FALSE],
                             kinship = kin[[18]], addcovar = add_covar)

####################
# Plot effects
quartz()
par(mfrow = c(2,1))
plot(x = qtl.delta_6.blup, map = map[[6]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 6", scan1_output = qtl.delta)
plot(x = qtl.delta_ptprz1_6.blup, map = map[[6]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 6 with Ptprz1 Covariate", scan1_output = qtl.delta_ptprz1)

quartz()
plot(x = qtl.delta_18.blup, map = map[[18]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 18", scan1_output = qtl.delta)
plot(x = qtl.delta_armc4_18.blup, map = map[[18]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 18", scan1_output = qtl.delta)
plot(x = qtl.delta_zfp438_18.blup, map = map[[18]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 18", scan1_output = qtl.delta)
plot(x = qtl.delta_zfp438_armc4_18.blup, map = map[[18]], columns = 1:8, col = CCcolors,
     main = "Delta Eigengene 18", scan1_output = qtl.delta)

quartz()
plot(x = qtl.ptprz1.blup, map = map[[6]], columns = 1:8, col = CCcolors, 
     main = "Ptprz1 6", scan1_output = qtl.ptprz1)
plot(x = qtl.armc4.blup, map = map[[18]], columns = 1:8, col = CCcolors, 
     main = "Armc4 18", scan1_output = qtl.ptprz1)
plot(x = qtl.zfp438.blup, map = map[[18]], columns = 1:8, col = CCcolors, 
     main = "Zfp438 18", scan1_output = qtl.ptprz1)

##############
# SNP Association Chromosome 6
##############
assoc.delta_ptprz1_6 <- assoc_mapping(probs = probs, pheno = phenos_genes, idx = 6, 
                                  addcovar = add_covar_ptprz1, k = kin, markers = snps, chr = 6,
                                  start = 15, end = 25, ncores = 4)

assoc.delta_6 <- assoc_mapping(probs = probs, pheno = phenos_genes, idx = 6, 
                                      addcovar = add_covar, k = kin, markers = snps, chr = 6,
                                      start = .01, end = 16, ncores = 4)

assoc.delta_18 <- assoc_mapping(probs = probs, pheno = phenos_genes, idx = 6, addcovar = add_covar, k = kin,
                                markers = snps, chr = 18, start = .01, end = 16, ncores = 4)
# Association mapping doesnt work anymore?!?!
# create gene matrix so I can plot genes under the snps
genes <- data.frame(chr = annot.mrna$chr, 
                    start = annot.mrna$start,
                    stop = annot.mrna$end, 
                    strand = annot.mrna$strand,
                    Name = annot.mrna$symbol, 
                    stringsAsFactors = FALSE)

# Plot snp associations
quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.delta_ptprz1_6[[1]], snpinfo = assoc.delta_ptprz1_6[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 6 & genes$start > 15e6 & genes$stop < 25e6,], 
           xlim = c(15, 25))

quartz()
par(mfrow = c(2,1))
plot_snpasso(scan1output = assoc.delta_6[[1]], snpinfo = assoc.delta_6[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 6 & genes$start > .01e6 & genes$stop < 16e6,], 
           xlim = c(.01, 16))

##### GENES UNDER CHR 6 (Ptprz1 covar) PEAK (@ approx. 23 Mbp)
# correlations: 
# Ing3, Gm26719, Cped1, Fam3c, Gm8927, Aass, Gm5989, Cadps2
ing3 <- get_exp_dat("Ing3")
gm26719 <- get_exp_dat("Gm26719")
cped1 <- get_exp_dat("Cped1")
fam3c <- get_exp_dat("Fam3c")
gm8927 <- get_exp_dat("Gm8927")
aass <- get_exp_dat("Aass")
gm5989 <- get_exp_dat("Gm5989")
cadps2 <- get_exp_dat("Cadps2")
wnt16 <- get_exp_dat("Wnt16")

cor(delta_eigengene, ing3)
# No
cor(delta_eigengene, gm26719)
# meh
cor(delta_eigengene, cped1)
# meh
cor(delta_eigengene, fam3c)
# No
cor(delta_eigengene, gm8927)
# No
cor(delta_eigengene, aass) # one good one
# -0.431
cor(delta_eigengene, gm5989)
# no
cor(delta_eigengene, cadps2)
# meh 
cor(delta_eigengene, wnt16)
# meh 
cor(delta_eigengene, ptprz1_exp)
# 0.8286


add_covar_ptprz1_aass <- model.matrix(~Sex + Generation + diet_days + ptprz1_exp + aass, data = temp)[,-1]

qtl.delta_ptprz1_aass <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                         kinship = kin, addcovar = add_covar_ptprz1_aass, cores = 4)

quartz()
plot_scan1(qtl.delta_ptprz1_aass, map = map, main = "Ptprz1 and Ptprz1, Aass", col = "red")
plot_scan1(qtl.delta_ptprz1, map = map, add = TRUE, col = "black")
# Aass DOES NOT MAKE THE PEAK DROP!!!

#########################################
# Genes under chr 6 peak (original at ~4Mbp)
tfpi2 <- get_exp_dat("Tfpi2")
gm26696 <- get_exp_dat("Gm26696")
casd1 <- get_exp_dat("Casd1")
ppp1r9a <- get_exp_dat("Ppp1r9a")
pon3 <- get_exp_dat("Pon3")
dync1i1 <- get_exp_dat("Dync1i1")

cor(tfpi2, delta_eigengene)
# - 0.08370911
cor(gm26696, delta_eigengene)
cor(casd1, delta_eigengene)
cor(ppp1r9a, delta_eigengene)
# HIGHEST CORELATION  0.2845801
cor(pon3, delta_eigengene)
cor(dync1i1, delta_eigengene)

add_covar_ppp1r9a <- model.matrix(~Sex + Generation + diet_days + ppp1r9a, data = temp)[,-1]

qtl.delta_ppp1r9a <- scan1(genoprobs = probs, pheno = phenos_genes[, colnames(phenos_genes) == "delta_eigengene", drop = FALSE],
                           kinship = kin, addcovar = add_covar_ppp1r9a, cores = 4)

quartz()
plot_scan1(qtl.delta_ppp1r9a, map = map, main = "Delta with ppp1r9a", col = "red")
plot_scan1(qtl.delta, map = map, add = TRUE, col = "black")

# I CANT FIND ANYTHING THAT DROPS THE LOD PEAK!! >:(

###################################
# REGRESSION OF Armc4, Zfp438
sex <- annot.samples$Sex
QD18 <- get_genoprob(18, 5.990333)

# Both genes and the peak have a significant effect on the eigengene
anova(lm(delta_eigengene~sex + armc4_exp)) # 6.884e-06 ***
anova(lm(delta_eigengene~sex + zfp438_exp)) # 0.0003357 ***
anova(lm(delta_eigengene~sex + QD18)) # 1.097e-06 ***

# Neither gene completely ablates the effect of the other
anova(lm(delta_eigengene~sex + armc4_exp + zfp438_exp)) 
# armc4_exp 6.274e-06 ***
# zfp438_exp 0.03637 * 
anova(lm(delta_eigengene~sex + zfp438_exp + armc4_exp))
# zfp438 0.0002759 ***
# armc4_exp 0.0006165 ***

# Neither gene completely ablates the effect of the peak genoprob
anova(lm(delta_eigengene~sex + armc4_exp + QD18)) 
# armc4_exp 4.767e-06 ***
# QD18 0.005129 **
anova(lm(delta_eigengene~sex + zfp438_exp + QD18)) 
# zfp438_exp 0.0002045 ***
# QD18 2.183e-05 ***

# Both genes together don't completely ablate the effect of the peak genoprob
anova(lm(delta_eigengene~sex + armc4_exp + zfp438_exp + QD18)) 
# armc4_exp 4.089e-06 *** 
# zfp438_exp 0.032693 * 
# QD18 0.002382 **  
anova(lm(delta_eigengene~sex + zfp438_exp + armc4_exp + QD18)) 
# zfp438_exp 0.0002070 ***
# armc4_exp 0.0004765 ***
# QD18 0.0023819 ** 



#################
# REGRESSION CHROMOSOME 6 PEAK
QD6 <- get_genoprob(6, 4.776124)
QD6_ptprz1 <- get_genoprob(6, 23.01854)

anova(lm(delta_eigengene~sex + QD6))
# QD6: 1.223e-07

# The genoprobs at the new peak dont even have a significant effect on the delta eigengene
anova(lm(delta_eigengene~sex + QD6_ptprz1))
# QD6_ptprz1: 0.1091 

anova(lm(delta_eigengene~sex + ptprz1_exp))
# ptprz1_exp: < 2.2e-16 ***

anova(lm(ptprz1_exp~sex + QD6))
# QD6: 6.213e-11 ***

anova(lm(ptprz1_exp~sex + QD6_ptprz1))
# QD6_ptprz1: 2.274e-11 ***

anova(lm(delta_eigengene~sex + ptprz1_exp + QD6))
# ptprz1_exp: <2e-16 ***
# QD6: 0.6164

anova(lm(delta_eigengene~sex + QD6 + ptprz1_exp))
# ptprz1_exp: <2e-16 ***
# QD6: <2e-16 ***
# Ptprz1_exp passes mediation analysis!

# Plotting LMs
plot(lm(delta_eigengene~sex))
plot(lm(delta_eigengene~weight_16wk))
plot(lm(delta_eigengene~weight_16wk + ins2_exp))

# Plot ptprz1_exp vs delta eigengene
quartz()
ggplot(phenos_genes, aes(x = delta_eigengene, y = ptprz1_exp, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
cor(ptprz1_exp, delta_eigengene)
# 0.8286802

#########################################
# Look at peak drops (you don't need to look at this anymore because you got intermediate to work)
max(qtl.delta, map = map, chr = 18)
# 6_4776124

qtl.delta["18_5990333",]
# 8.463
qtl.delta_armc4["18_5990333",]
# 3.029
qtl.delta_zfp438["18_5990333", ]
# 5.80415
qtl.delta_zfp438_armc4["18_5990333", ]
# 2.832454
#########################################

# Curious: do alpha cells mediate weight or does weight mediate alpha cells?
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow
QA11 <- get_genoprob(11, 16.65394)

# Is body weight driving the expression of the alpha eigengene through the chromosome 11 peak?
#i) alpha_eigengene is linked to alpha eigengene expression
anova(lm(alpha_eigengene~sex + weight_16wk)) # < 2.2e-16 ***
#ii) alpha_eigengene is linked to chr11 peak 
anova(lm(alpha_eigengene~sex + QA11)) #0.0002085 ***
#iii) alpha_eigengene not linked to QA11 after accounting for weight_16wk
anova(lm(alpha_eigengene~sex +  weight_16wk + QA11 )) 
# still linked (0.03781)
#iv) weight_16wk still linked to alpha eigengene after accounting for QA11 
anova(lm(alpha_eigengene~sex + QA11 + weight_16wk)) # Both significant

# OTHER WAY: Alpha cell is the driver
anova(lm(weight_16wk~sex + alpha_eigengene)) # < 2.2e-16 ***
anova(lm(weight_16wk~sex + QA11)) #0.005431 ** 
anova(lm(weight_16wk~sex +  alpha_eigengene + QA11 )) 
# alpha: <2e-16
# QA11: 0.4478
anova(lm(weight_16wk~sex + QA11 + alpha_eigengene)) # Both significant
# WORKS!!!!

#########################################
# One last try: look for candidate genes under the delta cell qtl peaks

# first: use intermediate to look for genes that significantly drop the peaks

#prepare to run intermediate
source("Intermediate_Scripts/plot.mediation.R")
source("Intermediate_Scripts/mediation.scan.R")
source("Intermediate_Scripts/gmb.coordinates.R")
load("data/mouse.chrlen.rda")

# Run mediation on chr 6 peak
med6 <- mediation.scan(target = delta_eigengene, mediator = rankz.mrna, annotation = annot.mrna, 
                       covar = add_covar, qtl.geno = QD6)
# Plot genes
quartz()
plot(med6)
# Klf14, Vwde, Mest, Kcnd2????

# Run mediation on chr 18 peak
med18 <- mediation.scan(target = delta_eigengene, mediator = rankz.mrna, annotation = annot.mrna, 
                        covar = add_covar, qtl.geno = QD18)
# Plot genes
quartz()
plot(med18)
# Armc4, Osbpl1a, Tmem241, Mib1


ptprz1_6 <- mediation.scan(target = ptprz1_exp, mediator = rankz.mrna, annotation = annot.mrna, 
                       covar = add_covar, qtl.geno = QD6)
# Plot genes
quartz()
plot(ptprz1_6)
# no mediators besides itself

#######################
# CURIOUS
med_bw_11 <- mediation.scan(target = weight_16wk, mediator = rankz.mrna, annotation = annot.mrna, 
                        covar = add_covar, qtl.geno = QA11)
quartz()
plot(med_bw_11)

med_alpha_11 <- mediation.scan(target = alpha_eigengene, mediator = rankz.mrna, annotation = annot.mrna, 
                            covar = add_covar, qtl.geno = QA11)
quartz()
plot(med_alpha_11)
####################

# I guess I just end with the fact that I found some candidate genes...






















