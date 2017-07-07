# Result Confirmation
# Sadie Allen
# July 7, 2017
# In this script I will rerun several analyses to confirm results I have obtained

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

##### Confirm phenotype data is normalized #####
colnames(ghrelin_list)

hist(ghrelin_list$food_ave)
hist(log10(ghrelin_list$food_ave)) # This looks more normal... I should save it as this... right?

hist(ghrelin_list$weight_sac)
hist(log10(ghrelin_list$weight_sac)) # This one also looks more normal after being log transformed

hist(ghrelin_list$weight_6wk)
hist(log10(ghrelin_list$weight_6wk)) # This one is also more normal

hist(ghrelin_list$weight_change)
hist(log10(ghrelin_list$weight_change)) # Neither of these seem super normal so I will leave it how it is

# Normalize phenotypes
ghrelin_list$food_ave <- log10(ghrelin_list$food_ave)
ghrelin_list$weight_sac <- log10(ghrelin_list$weight_sac)
ghrelin_list$weight_6wk <- log10(ghrelin_list$weight_6wk)

ghrelin_list$X <- NULL

# Save this to be used in subsequent analyses
write.csv(ghrelin_list, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")

ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")

# Rerun important QTL scans
# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ghsr, map = map, threshold = 6, drop = 1.5)
# Peak: chr 18, pos 4.474313, ci: 0.026149-10.21129
Q18 <- get_genoprob(18, 4.474313)

# call important gene expressions
svil_exp <- ghrelin_list$svil_exp
fzd8_exp <- ghrelin_list$fzd8_exp
mpp7_exp <- ghrelin_list$mpp7_exp
ccny_exp <- ghrelin_list$ccny_exp
mtpap_exp <- ghrelin_list$mtpap_exp
efnb3_exp <- ghrelin_list$efnb3_exp

### CONFIRMING RESULTS OF MEDIATION ANALYSIS: 

# Food Intake as a Mediator of the Relationship Between Ghsr Expression and Body Weight
food_ave <- ghrelin_list$food_ave
ghsr_exp <- ghrelin_list$ghsr_exp
weight_sac <- ghrelin_list$weight_sac
sex <- ghrelin_list$sex
# i) weight_sac is linked to ghsr_exp
anova(lm(weight_sac ~ sex + ghsr_exp))
# p = 9.121e-15 CONDITION MET
# ii) food_ave is linked to ghsr_exp
anova(lm(food_ave ~ sex + ghsr_exp))
# p = 1.313e-13 CONDITION MET
# iii) weight_sac not linked to ghsr_exp after accounting for food_ave
anova(lm(weight_sac ~ sex + food_ave + ghsr_exp))
# p = 4.749e-05 CONDITION NOT MET
# iv) food_ave still linked to ghsr_exp after accounting for weight_sac
anova(lm(food_ave ~ sex + weight_sac + ghsr_exp))
# p = 7.69e-04 CONDITION MET
# As before, food_ave did not meet the conditions of being a mediator 

# Ghsr expression as a mediator of the relationship between Q18 genotype and food intake
# i) food_ave is linked to Q18
anova(lm(food_ave ~ sex + Q18))
# p = 0.6038 CONDITION NOT MET
# ii) ghsr_exp is linked to Q18
anova(lm(ghsr_exp ~ sex + Q18))
# p = 8.392e-07 CONDITION MET
# iii) food_ave not linked to Q18 after accounting for ghsr_exp
anova(lm(food_ave ~ sex + ghsr_exp + Q18))
# p = 0.3912 CONDITION MET
# iv) ghsr_exp still linked to Q18 after accounting for food_ave
anova(lm(ghsr_exp ~ sex + food_ave + Q18))
# p = 3.755e-07 CONDITION MET
# Since Q18 is not linked to food_ave, there can be no mediator between the two

# Svil expression as a mediator of the relationship between Q18 genotype and ghsr_exp
# i) ghsr_exp is linked to Q18
anova(lm(ghsr_exp ~ sex + Q18))
# P = 8.392e-07 CONDITION MET
# ii) svil_exp is linked to Q18
anova(lm(svil_exp ~ sex + Q18))
# P < 2.2e-16 CONDITION MET
# iii) ghsr_exp not linked to Q18 after accounting for svil_exp
anova(lm(ghsr_exp ~ sex + svil_exp + Q18))
# P = 0.0003026 CONDITION NOT MET
# iv) svil_exp still linked to Q18 after accounting for ghsr_exp
anova(lm(svil_exp ~ sex + ghsr_exp + Q18))
# P < 2.2e-16
# Svil did not meet the third condition of mediation analysis, but it seemed close

# Fzd8 expression as a mediator of the relationship between Q18 genotype and ghsr_exp
# i) ghsr_exp is linked to Q18
anova(lm(ghsr_exp ~ sex + Q18))
# P = 8.392e-07 CONDITION MET
# ii) fzd8_exp is linked to Q18
anova(lm(fzd8_exp ~ sex + Q18))
# P = 3.306e-07 CONDITION MET
# iii) ghsr_exp not linked to Q18 after accounting for fzd8_exp
anova(lm(ghsr_exp ~ sex + fzd8_exp + Q18))
# P = 3.354e-07 CONDITION NOT MET
# iv) fzd8_exp still linked to Q18 after accounting for ghsr_exp
anova(lm(fzd8_exp ~ sex + ghsr_exp + Q18))
# P = 1.318e-07 CONDITION MET
# Fzd8 did not meet the third condition of mediation analysis

# Efnb3 expression as a mediatior of the relationship between Ghsr expression and food intake
# i) food_ave is linked to ghsr_exp
anova(lm(food_ave ~ sex + ghsr_exp))
# P = 1.313e-13
# Efnb3 expression is linked to Ghsr expression
anova(lm(efnb3_exp ~ sex + ghsr_exp))
# P < 2.2e-16
# food intake NOT linked to Ghsr expression after accounting for Efnb3
anova(lm(food_ave ~ sex + efnb3_exp + ghsr_exp)) # a weaker link
# P = 0.03282
# Efnb3 expression and Ghsr expression linked after accounting for food intake
anova(lm(efnb3_exp ~ sex + food_ave + ghsr_exp))
# P = 2.2e-16
# Efnb3 did not meet the third condition of mediation analysis, but it was close

### CONFIRMING RESULTS OF BIC ANALYSIS

# Determining the nature of the relationship between food_ave, ghsr_exp, and Q18
triple.fit(ghsr_exp, food_ave, Q18)
# independent: 175.44958
# reactive: 71.01257
# causal: 48.01773 --> BIC analysis suggests a causal model! Q18 -> ghsr_exp -> food_intake
# complex: 76.62758

# BIC Modeling for genes under chromosome 18 peak 
# Here I made a mistake earlier because I used the genotype at the peak for the specific gene 
# instead of using the ghsr distant qtl peak (although these should be similar, they are not identical)
# genes: Svil, Fzd8, Mpp7, Ccny, Mtpap
# In each of these triple fits I will be testing the relationship between ghsr expression, GOI expression,
# and the chromosome 18 peak
## Svil ##
triple.fit(svil_exp, ghsr_exp, Q18)
# independent: 1994.637
# reactive: 2146.078
# causal: 1980.042 --> BIC analysis suggests a causal model! Q18 -> svil_exp -> ghsr_exp
# complex: 2000.137
## Fzd8 ##
triple.fit(fzd8_exp, ghsr_exp, Q18)
# independent: 2142.582
# reactive: 2129.404
# causal: 2113.792 --> BIC analysis suggests a causal model! Q18 -> fzd8 -> ghsr_exp
# complex: 2125.054
## Mpp7 ##
triple.fit(mpp7_exp, ghsr_exp, Q18)
# independent: 2119.477 
# reactive: 2115.104
# causal: 2076.306
# complex: 2072.883 --> BIC analysis suggests a complex model
## Ccny ##
triple.fit(ccny_exp, ghsr_exp, Q18)
# independent: 2156.927 
# reactive: 2116.839
# causal: 2115.092
# complex: 2139.345 --> BIC analysis does not point strongly to any one model
triple.fit(mtpap_exp, ghsr_exp, Q18)
# independent: 2157.644 
# reactive: 2127.124
# causal: 2122.569
# complex: 2143.367 --> BIC analysis does not point strongly to any one model
# So, strong evidence for causal models for Fzd8 and Svil, but not for any of the other genes

### SIMPLE BIC MODELING ###
# Svil expression affects ghsr expression
BIC(lm(ghsr_exp~sex))
# 1009.857
BIC(lm(ghsr_exp~ sex + svil_exp))
#997.6577
#Q18 affects svil expression
BIC(lm(svil_exp ~ sex))
#1065.1
BIC(lm(svil_exp ~ sex + Q18))
#916.3928
#Q18 affects ghsr expression
# WHY?!?!?!
BIC(lm(ghsr_exp~sex))
# 1009.857
BIC(lm(ghsr_exp~sex + Q18))
#1009.763
BIC(lm(ghsr_exp~sex + svil_exp + Q18))
#1011.416
# This seemed weird the first time, and now I have run this code twice with the same result and
# I have no idea what I am doing wrong...

# food_ave is affected by ghsr_exp and efnb3_exp
BIC(lm(food_ave~sex))
# 472.1798 # NOW -1016
BIC(lm(food_ave~sex + ghsr_exp))
# 423.0559 # -1065
BIC(lm(food_ave~sex + ghsr_exp + efnb3_exp))
# 409.1854 # -1079










