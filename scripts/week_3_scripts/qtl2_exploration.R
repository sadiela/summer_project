# QTL2 Exploration
# Sadie Allen
# June 22, 2017
# In this script, I work through some of the functions in the
# R/qtl2 package, created by Karl Broman, to gain a better 
# understanding of the analyses I will be performing. My script
# will be based off of Broman's qtl2 user guide, Daniel Gatti's
# qtl2 demo, and a QTL mapping script created by Petr Simecek.

# Load necessary libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)

# The first step in the user guide is to calculate genotype 
# probabilities using the calc_genoprob() function. However,
# I do not have to do this because when I load the islet data:
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# It includes an object called "genoprobs" that contains
# this information. However, I do need to convert it from DOQTL
# to qtl2 format:
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# (While I'm at it, I might as well add some phenotypical data
# that I will need later)
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
# Remove some phenotypes from my dataframe
mod_ghrelin_shortlist <- ghrelin_shortlist
mod_ghrelin_shortlist$sex <- NULL
mod_ghrelin_shortlist$G83_ins_secrete <- NULL
mod_ghrelin_shortlist$G167_ins_secrete <- NULL
mod_ghrelin_shortlist$X <- NULL
mod_ghrelin_shortlist$ghrs <- NULL
mod_ghrelin_shortlist$sstr <- NULL
mod_ghrelin_shortlist$ins2 <- NULL
mod_ghrelin_shortlist$ghrl <- NULL
mod_ghrelin_shortlist$ins1 <- NULL
mod_ghrelin_shortlist$sst <- NULL
mod_ghrelin_shortlist$Mouse.ID <- NULL

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"


## Calculate a Kinship Matrix ##
# The next suggestion in the user guide is to calculate a kinship
# matrix, which is necessary if you want to perform a genome scan
# by a linear mixed model. It accounts for the relationships among
# individuals. Kinship can be calculated using the calc_kinship()
# function and takes genotype probabilities as input.  
kin <- calc_kinship(probs = probs, type = "loco", cores = 4)

## Special Covariates for the X Chromosome ##
# The user guide then discusses special covariates that may need to 
# be included under the null hypothesis of no QTL to avoid spurious
# evidence of linkage. Broman uses the get_x_covar() function to do 
# this. Petr's script does not do this step. In addition, I am not
# sure what the argument it takes is, so I am not going to do this
# step either. 

## Additive Covariates ##
# Both Dan and Petr create an additive covariate object that 
# incorporates sex, generation, and diet days. Since I want to
# include diet days as a covariate but it is not in the annot.samples
# object, I will add it from my phenotype dataframe:
#Mouse.ID <- as.character(ghrelin_shortlist$Mouse.ID)
#diet_days_id <- data.frame(Mouse.ID = Mouse.ID, diet_days =ghrelin_shortlist$diet_days)
temp = merge(annot.samples, mod_ghrelin_shortlist, by = "row.names")
rownames(temp) = temp[,1]
# annot.samples <- merge(annot.samples, diet_days_id, by.x="Mouse.ID", by.y = "Mouse.ID")
# Now, create the covariate
add_covar <- model.matrix(~Sex + Generation + diet_days, data = temp)[,-1]
# Question: could I also use the covar object here? It contains the same
# information, and if I'm not using it here, when would I?

# remove diet days from pheno list
mod_ghrelin_shortlist$diet_days <- NULL

## Miscellaneous Steps ##
#convert a marker map organized as data frame to a list
map <- map_df_to_list(map = snps, pos_column = "bp")
#make sure the rownames of the expression data and phenotype data match
#rownames(mod_ghrelin_shortlist) <- rownames(expr.mrna)

## Performing a Genome Scan ##
# To perform a genome scan by Haley-Knott regression, use the function
# scan1() in qtl2scan. The inputs of this function are the genotype
# probabilities, a matrix of phenotypes, optional additive and interactive
# covariates, and the X chromosome covariates. This computer (MLG-CCMBP01)
# does not have a lot of computational power, so I will run a scan of
# one phenotype at a time. 
qtl.glu0min <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,1, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.glusac <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,2, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.ins0min <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,3, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.inssac <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,4, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.foodave <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,5, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.weightsac <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,6, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.g33inssecrete <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,7, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.numislets <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,8, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.insperislet <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,9, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.wpic <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,10, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

qtl.weightchangeave <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,11, drop = FALSE],
                     kinship = kin, addcovar = add_covar, cores = 4)

## Plotting QTL Scans ##
#convert a marker map organized as data frame to a list
plot_scan1(x = qtl.glu0min, map = map, main = colnames(qtl.glu0min)[1])
plot_scan1(x = qtl.glusac, map = map, main = colnames(qtl.glusac)[1])
plot_scan1(x = qtl.ins0min, map = map, main = colnames(qtl.ins0min)[1])
plot_scan1(x = qtl.inssac, map = map, main = colnames(qtl.inssac)[1])
plot_scan1(x = qtl.foodave, map = map, main = colnames(qtl.foodave)[1])
plot_scan1(x = qtl.g33inssecrete, map = map, main = colnames(qtl.g33inssecrete)[1])
plot_scan1(x = qtl.weightsac, map = map, main = colnames(qtl.weightsac)[1])
plot_scan1(x = qtl.wpic, map = map, main = colnames(qtl.wpic)[1])
plot_scan1(x = qtl.numislets, map = map, main = colnames(qtl.numislets)[1])
plot_scan1(x = qtl.insperislet, map = map, main = colnames(qtl.insperislet)[1])
# Insperislet has an interesting peak on chr 5!
plot_scan1(x = qtl.weightchangeave, map = map, main = colnames(qtl.weightchangeave)[1])
# How can I put labels on my qtl scans??

## Find Peaks ##
find_peaks(qtl.glu0min, map, threshold = 5, drop = 1.5)
## chr 7, pos 78.04246
## chr 8, pos 109.53512 (highest peak)
## chr 12, pos 34.66649
find_peaks(qtl.glusac, map, threshold = 5, drop = 1.5)
## chr 3, 152
## chr 13, 71 (highest)
## chr 15, 100
## chr 18, 68
find_peaks(qtl.ins0min, map, threshold = 6, drop =1.5)
## 11, 81 (highest)
## 12, 77
## 13, 97
## 17, 73
find_peaks(qtl.inssac, map, threshold = 6, drop =1.5)
## 11, 73 (highest)
## 17, 69
## 18, 75
## 20, 47
find_peaks(qtl.foodave, map, threshold = 6, drop =1.5)
## 1, 172
## 7, 136 (highest)
find_peaks(qtl.g33inssecrete, map, threshold = 5, drop =1.5)
## 1, 51
## 3, 15
## 5, 125
## 14, 92 (highest)
find_peaks(qtl.weightsac, map, threshold = 5.5, drop =1.5)
## 3, 151
## 9, 103 (highest)
## 17, 31
## 20, 48
find_peaks(qtl.wpic, map, threshold = 5.5, drop =1.5)
## 8, 22
find_peaks(qtl.numislets, map, threshold = 6, drop =1.5)
## 18, 58
find_peaks(qtl.insperislet, map, threshold = 6, drop =1.5)
## 5, 91 (highest)
## 17, 21 
find_peaks(qtl.weightchangeave, map, threshold = 6, drop =1.5)
## 17, 31 (highest)
## 20, 139

## Calculate Founder Allele Effects (Coefficients) and Blups ##
# Glu_0min 
chr = 8
qtl.glu0min.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,1, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.glu0min.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.glu0min)[1], scan1_output = qtl.glu0min)

qtl.glu0min.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,1, drop = FALSE],
                                 kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.glu0min.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.glu0min)[1], scan1_output = qtl.glu0min)


# Glu_sac
chr = 13
qtl.glusac.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,2, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.glusac.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.glusac)[1], scan1_output = qtl.glusac)

qtl.glusac.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,2, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.glusac.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.glusac)[1], scan1_output = qtl.glusac)

# Ins_0min
chr = 11
qtl.ins0min.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,3, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.ins0min.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ins0min)[1], scan1_output = qtl.ins0min)

qtl.ins0min.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,3, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.ins0min.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.ins0min)[1], scan1_output = qtl.ins0min)

# Ins_sac
chr = 11
qtl.inssac.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,4, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.inssac.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.inssac)[1], scan1_output = qtl.inssac)

qtl.inssac.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,4, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.inssac.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.inssac)[1], scan1_output = qtl.inssac)

# Food_ave
chr= 7
qtl.foodave.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,5, drop = FALSE],
                             kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.foodave.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.foodave)[1], scan1_output = qtl.foodave)

qtl.foodave.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,5, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.foodave.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.foodave)[1], scan1_output = qtl.foodave)

# G33_ins_secrete 
chr= 14
qtl.g33inssecrete.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,6, drop = FALSE],
                              kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.g33inssecrete.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.g33inssecrete)[1], scan1_output = qtl.g33inssecrete)

qtl.g33inssecrete.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,6, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.g33inssecrete.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.g33inssecrete)[1], scan1_output = qtl.g33inssecrete)

# Weight_sac
chr= 9
qtl.weightsac.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,7, drop = FALSE],
                                    kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.weightsac.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weightsac)[1], scan1_output = qtl.weightsac)

qtl.weightsac.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,7, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.weightsac.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weightsac)[1], scan1_output = qtl.weightsac)

# Num_islets
chr= 18
qtl.numislets.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,8, drop = FALSE],
                                kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.numislets.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.numislets)[1], scan1_output = qtl.numislets)

qtl.numislets.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,8, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.numislets.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.numislets)[1], scan1_output = qtl.numislets)

# Ins_per_islet
chr = 5
qtl.insperislet.coef = scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,9, drop = FALSE],
                                 kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.insperislet.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.insperislet)[1], scan1_output = qtl.insperislet)

qtl.insperislet.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,9, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.insperislet.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.insperislet)[1], scan1_output = qtl.insperislet)

# WPIC
chr= 8
qtl.wpic.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,10, drop = FALSE],
                                kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.wpic.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.wpic)[1], scan1_output = qtl.wpic)

qtl.wpic.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,10, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.wpic.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.wpic)[1], scan1_output = qtl.wpic)

# weight_change_ave
chr= 17
qtl.weightchangeave.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,11, drop = FALSE],
                                kinship = kin[[chr]], addcovar = add_covar)
plot(x = qtl.weightchangeave.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weightchangeave)[1], scan1_output = qtl.weightchangeave)

qtl.weightchangeave.blup = scan1blup(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,11, drop = FALSE],
                             kinship=kin[[chr]], addcovar = add_covar)
plot(x = qtl.weightchangeave.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl.weightchangeave)[1], scan1_output = qtl.weightchangeave)

## Association Mapping? ##







