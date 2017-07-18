# Association Mapping Verification
# Sadie Allen
# July 13, 2017
# Verifying mediation results and looking for snips

rm(list = ls())

#Load Libraries
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

arg1_exp <- get_exp_dat("Arg1")
ptprz1_exp <- get_exp_dat("Ptprz1")
kcnd2_exp <- get_exp_dat("Kcnd2")
cacna1h_exp <- get_exp_dat("Cacna1h")
armc4_exp <- get_exp_dat("Armc4")
dscam_exp <- get_exp_dat("Dscam")
spock3_exp <- get_exp_dat("Spock3")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
hhex_exp <- get_exp_dat("Hhex")
sst_exp <- get_exp_dat("Sst")
slc16a7_exp <- get_exp_dat("Slc16a7")
efnb3_exp <- get_exp_dat("Efnb3")
ly6h_exp <- get_exp_dat("Ly6h")
wnt10a_exp <- get_exp_dat("Wnt10a")
f5_exp <- get_exp_dat("F5")
ncam2_exp <- get_exp_dat("Ncam2")
rbm12b2_exp <- get_exp_dat("Rbm12b2")

ghrelin_assoc <- cbind(ghrelin_list, arg1_exp, ptprz1_exp, kcnd2_exp, cacna1h_exp, armc4_exp,
                       sst_exp, dscam_exp, slc16a7_exp, rbm12b2_exp)

# Prepare to run QTL scans
load("data/qtl_prep.RData")
# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# QTLS I need:
# Chr 18: svil, armc4, cacna1h, ghsr, efnb3
# Chr 6: cacna1h, slc16a7, ptprz1
# chr 4: arg1, ghsr, efnb3

qtl.arg1 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,16, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.arg1, map = map, threshold = 6, drop = 1.5)
#lodindex lodcolumn chr      pos      lod    ci_lo    ci_hi
#1        1  arg1_exp  10 27.55391 8.598977 25.64208 31.20788
#2        1  arg1_exp  16 91.45495 6.624899 91.32522 93.08095

qtl.ptprz1 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,17, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ptprz1, map = map, threshold = 6, drop = 1.5)
#1        1 ptprz1_exp   1 92.45837  7.082260 91.94546 93.15224
#2        1 ptprz1_exp   4 12.77392  6.325366 11.68454 85.18468
#3        1 ptprz1_exp   6 23.20927 17.716343 21.38721 23.34616
#4        1 ptprz1_exp   9 65.58070  8.093590 65.42280 67.51487

qtl.slc16a7 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,23, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.slc16a7, map = map, threshold = 6, drop = 1.5)
#lodindex   lodcolumn chr        pos       lod      ci_lo     ci_hi
#1        1 slc16a7_exp   6   4.776124  9.123031   0.738935  13.22898
#2        1 slc16a7_exp  10 125.290314 11.902760 124.933768 127.08662
#3        1 slc16a7_exp  13  98.487897  6.016090  59.557325  99.03916

qtl.efnb3 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,15, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.efnb3, map = map, threshold = 6, drop = 1.5)
#lodindex lodcolumn chr       pos      lod     ci_lo     ci_hi
#1        1 efnb3_exp   4  5.566854 7.127355  0.207518  13.32408
#2        1 efnb3_exp  11 19.230549 6.541782 17.680768  19.61891
#3        1 efnb3_exp  15 71.061095 6.017569 59.966210 102.66890
#4        1 efnb3_exp  17 79.310467 6.262556 25.947374  79.53880
#5        1 efnb3_exp  18  3.351052 7.308181  0.026149  10.25686

qtl.svil <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,10, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.svil, map = map, threshold = 6, drop = 1.5)
#lodindex lodcolumn chr       pos       lod     ci_lo     ci_hi
#1        1  svil_exp   9 75.921080  6.249971 74.813305 80.130566
#2        1  svil_exp  18  4.772599 49.236011  4.431701  4.857823

qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,7, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ghsr, map = map, threshold = 6, drop = 1.5)
#  lodindex lodcolumn chr       pos      lod    ci_lo    ci_hi
#1        1  ghsr_exp   4 11.984686 6.859136 0.207518 13.32408
#2        1  ghsr_exp  18  4.474313 7.822138 0.026149 10.21129
qtl.armc4 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,20, drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.armc4, map = map, threshold = 6, drop = 1.5)
#lodindex lodcolumn chr        pos       lod      ci_lo      ci_hi
#1        1 armc4_exp   4 124.519326   6.25769 119.964458 125.605186
#2        1 armc4_exp  18   7.761497 120.69338   7.003513   8.402237

qtl.cacna1h <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,19, drop = FALSE],
                           kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.cacna1h, map = map, threshold = 6, drop = 1.5)
#lodindex   lodcolumn chr        pos       lod      ci_lo      ci_hi
#1        1 cacna1h_exp   4   4.570163  6.389070   0.207518  13.255308
#2        1 cacna1h_exp   6   4.323241  8.036235   0.738935  11.331688
#3        1 cacna1h_exp  10 119.143146  6.672532 116.095387 119.944193
#4        1 cacna1h_exp  18   4.516925 10.112666   0.026149   8.530216
qtl.rbm12b2 <- scan1(genoprobs = probs, pheno = ghrelin_assoc[,24, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.rbm12b2, map = map, threshold = 6, drop = 1.5)
#lodindex   lodcolumn chr      pos      lod    ci_lo    ci_hi
#1        1 rbm12b2_exp   4 12.23352 16.32164 11.68454 12.63638

plot_scan1(x = qtl.ghsr, map = map, main = "ghsr", col = "black")
plot_scan1(x = qtl.efnb3, map = map, main = "efnb3", col = "black")
plot_scan1(x = qtl.cacna1h, map = map, main = "cacna1h", col = "black")
plot_scan1(x = qtl.rbm12b2, map = map, main = "rmb12b2")

plot_scan1(x = qtl.ptprz1, map = map, main = "ptprz1")
plot_scan1(x = qtl.slc16a7, map = map, main = "slc16a7")



# Plotting
### CHROMOSOME 18
quartz()
par(mfrow = c(5, 1))
plot_scan1(x = qtl.ghsr, map = map, main = "Ghsr_exp")
plot_scan1(x = qtl.efnb3, map = map, main = "Efnb3_exp")
plot_scan1(x = qtl.cacna1h, map = map, main = "Cacna1h_exp")
plot_scan1(x = qtl.svil, map = map, main = "Svil_exp")
plot_scan1(x = qtl.armc4, map = map, main = "Armc4_exp")

# COEFFICIENT PLOTS (BLUPS)
chr = 18
qtl.ghsr.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,7, drop = FALSE],
                                 kinship = kin[[chr]], addcovar = add_covar)
qtl.efnb3.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,15, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
qtl.cacna1h.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,19, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
qtl.svil.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,10, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
qtl.armc4.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,20, drop = FALSE],
                           kinship = kin[[chr]], addcovar = add_covar)
quartz()
plot(x = qtl.ghsr.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Ghsr", scan1_output = qtl.ghsr)
plot(x = qtl.efnb3.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Efnb3", scan1_output = qtl.efnb3)
plot(x = qtl.cacna1h.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Cacna1h", scan1_output = qtl.cacna1h)
plot(x = qtl.svil.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Svil", scan1_output = qtl.svil)
plot(x = qtl.armc4.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Armc4", scan1_output = qtl.armc4)

assoc.ghsr <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 7, 
                            addcovar = add_covar, k = kin, markers = snps, chr = 18,
                            start = 0.207518, end = 13.32408, ncores = 4)
assoc.efnb3 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 15, 
                               addcovar = add_covar, k = kin, markers = snps, chr = 18,
                               start = 0.207518, end = 13.32408, ncores = 4)
assoc.armc4 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 20, 
                            addcovar = add_covar, k = kin, markers = snps, chr = 18,
                            start = 0.207518, end = 13.32408, ncores = 4)
assoc.cacna1h_18 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 19, 
                            addcovar = add_covar, k = kin, markers = snps, chr = 18,
                            start = 0.207518, end = 13.32408, ncores = 4)

genes <- data.frame(chr = annot.mrna$chr, 
                    start = annot.mrna$start,
                    stop = annot.mrna$end, 
                    strand = annot.mrna$strand,
                    Name = annot.mrna$symbol, 
                    stringsAsFactors = FALSE)

########################################
quartz()
par(mfrow = c(5,1))
plot_snpasso(scan1output = assoc.ghsr[[1]], snpinfo = assoc.ghsr[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.efnb3[[1]], snpinfo = assoc.efnb3[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.cacna1h[[1]], snpinfo = assoc.cacna1h[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.armc4[[1]], snpinfo = assoc.armc4[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 18 & genes$start > 0.207518e6 & genes$stop < 13.32408e6,], 
           xlim = c(0.207518, 13.32408))
########################################


quartz()
par(mfrow = c(3, 1))
plot_scan1(x = qtl.cacna1h, map = map, main = "Cacna1h_exp")
plot_scan1(x = qtl.slc16a7, map = map, main = "Slc16a7_exp")
plot_scan1(x = qtl.ptprz1, map = map, main = "Ptprz1_exp")

chr = 6
qtl.cacna1h_6.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,19, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
qtl.slc16a7.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,23, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
qtl.ptprz1.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,17, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)

quartz()
par(mfrow = c(3, 1))
plot(x = qtl.cacna1h_6.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Cacna1h", scan1_output = qtl.cacna1h)
plot(x = qtl.slc16a7.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Slc16a7", scan1_output = qtl.slc16a7)
plot(x = qtl.ptprz1.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Ptprz1", scan1_output = qtl.ptprz1)

assoc.cacna1h_6 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 19, 
                                 addcovar = add_covar, k = kin, markers = snps, chr = 6,
                                 start = 0.738935, end = 23.34616, ncores = 4)
assoc.slc16a7 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 23, 
                            addcovar = add_covar, k = kin, markers = snps, chr = 6,
                            start = 0.738935, end = 23.34616, ncores = 4)
assoc.ptprz1 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 17, 
                             addcovar = add_covar, k = kin, markers = snps, chr = 6,
                             start = 0.738935, end = 23.34616, ncores = 4)

quartz()
par(mfrow = c(5,1))
plot_snpasso(scan1output = assoc.cacna1h_6[[1]], snpinfo = assoc.cacna1h_6[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.slc16a7[[1]], snpinfo = assoc.slc16a7[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.ptprz1[[1]], snpinfo = assoc.ptprz1[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 6 & genes$start > 0.738935e6 & genes$stop < 23.34616e6,], 
           xlim = c(0.738935, 23.34616))

### CHROMOSOME 4
quartz()
par(mfrow = c(5, 1))
plot_scan1(x = qtl.ghsr, map = map, main = "Ghsr_exp")
plot_scan1(x = qtl.efnb3, map = map, main = "Efnb3_exp")
plot_scan1(x = qtl.cacna1h, map = map, main = "Cacna1h_exp")
plot_scan1(x = qtl.arg1, map = map, main = "Arg1_exp")
plot_scan1(x = qtl.rbm12b2, map = map, main = "Rbm12b2_exp")

chr = 4
qtl.ghsr_4.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,7, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)
qtl.efnb3_4.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,15, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)
qtl.cacna1h_4.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,19, drop = FALSE],
                                kinship = kin[[chr]], addcovar = add_covar)
qtl.arg1.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,16, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
qtl.rbm12b2.blup <- scan1blup(genoprobs = probs[,chr], pheno = ghrelin_assoc[,24, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)

quartz()
plot(x = qtl.ghsr_4.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Ghsr", scan1_output = qtl.ghsr)
plot(x = qtl.efnb3_4.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Efnb3", scan1_output = qtl.efnb3)
plot(x = qtl.cacna1h_4.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Cacna1h", scan1_output = qtl.cacna1h)
plot(x = qtl.arg1.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Arg1", scan1_output = qtl.arg1)
plot(x = qtl.rbm12b2.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = "Rbm12b2", scan1_output = qtl.rbm12b2)

assoc.ghsr_4 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 7, 
                               addcovar = add_covar, k = kin, markers = snps, chr = chr,
                               start = 0.207518, end = 13.32408, ncores = 4)
assoc.efnb3_4 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 15, 
                              addcovar = add_covar, k = kin, markers = snps, chr = chr,
                              start = 0.207518, end = 13.32408, ncores = 4)
assoc.cacna1h_4 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 19, 
                                 addcovar = add_covar, k = kin, markers = snps, chr = chr,
                                 start = 0.207518, end = 13.32408, ncores = 4)
assoc.arg1 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 16, 
                               addcovar = add_covar, k = kin, markers = snps, chr = chr,
                               start = 0.207518, end = 13.32408, ncores = 4)
assoc.rbm12b2 <- assoc_mapping(probs = probs, pheno = ghrelin_assoc, idx = 24, 
                              addcovar = add_covar, k = kin, markers = snps, chr = chr,
                              start = 0.207518, end = 13.32408, ncores = 4)


quartz()
par(mfrow = c(3,1))
plot_snpasso(scan1output = assoc.ghsr_4[[1]], snpinfo = assoc.ghsr_4[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.efnb3_4[[1]], snpinfo = assoc.efnb3_4[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.cacna1h_4[[1]], snpinfo = assoc.cacna1h_4[[2]], 
             drop.hilit = 1)
quartz
par(mfrow = c(3,1))
plot_snpasso(scan1output = assoc.arg1[[1]], snpinfo = assoc.arg1[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.rbm12b2[[1]], snpinfo = assoc.rbm12b2[[2]], 
             drop.hilit = 1)
plot_genes(genes[genes$chr == 4 & genes$start > 0.207518e6 & genes$stop < 13.32408e6,], 
           xlim = c(0.207518, 13.32408))

save.image(file = "data/qtls_blups_assoc.RData")

# Look at... stuff idk Dan showed me this
local.map  = qtl2plot:::snpinfo_to_map(assoc.cacna1h[[2]])
local.snps = qtl2plot:::expand_snp_results(assoc.cacna1h[[1]], local.map, assoc.cacna1h[[2]])






















