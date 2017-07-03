# Attempting linear models
# Sadie Allen
# June 27, 2017
# In this script I will be trying to adapt the triplefit function used on the BTBR 
# f2 cross for use in the DO dataset. I will also perform some mediation analysis

rm(list = ls())

## Load Libraries ##
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)

## Load Data ##
# Phenotype Data
ghrelin_list <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list[,1]
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load all functions
source("scripts/functions.R")

## Ghrelin & Food Intake ##
genotype_ghsr <- get_genoprob(18, 4.474313)
ghsr_exp <- ghrelin_list$ghsr_exp
food_ave <- ghrelin_list$food_ave
weight_sac <- ghrelin_list$weight_sac

triple.fit(ghsr_exp, food_ave, genotype_ghsr)
# independent: 175.06208
# reactive: 71.45350
# causal: 47.18930 Suggests a causal relationship! genotype -> gene expression -> phenotype
# complex: 77.67778

# MEDIATION ANALYSIS
#####
# check 4 conditions for Pdrg1 gene expression 
# as a mediator of Q2 effect on insulin
####  i) Insulin is linked to Q2
#anova(lm(INS.10wk ~ Sex + Q2, data = f2g$pheno))
# significant
####  ii) Pdrg1 gene expression is linked to Q2
#anova(lm(pdrg1_islet ~ Sex + Q2, data = f2g$pheno))
# significant
####  iii) Insulin not linked after accounting for Q2
#anova(lm(INS.10wk ~ Sex + pdrg1_islet + Q2, data = f2g$pheno))
# not significant * .
####  iv) Pdrg1 gene expression is still linked after 
# accounting for insulin
#anova(lm(pdrg1_islet ~ Sex + INS.10wk + Q2, data = f2g$pheno))
# significant ***
# all 4 conditions for a mediator are satisfied



# check conditions for ghsr_exp as a mediator of genotype_ghsr effect on food_ave
# i) food_ave linked to genotype_ghsr
anova(lm(food_ave~ghrelin_list$sex + genotype_ghsr))
# Not significant! >:(
# UNSUCCESSFUL

# check conditions for food_ave as a mediator of ghsr_exp effect on weight_sac
# i) weight_sac linked to ghsr_exp
anova(lm(weight_sac~ghrelin_list$sex + ghsr_exp))
# Significant! 
# ii) food_ave is linked to ghsr_exp
anova(lm(food_ave~ghrelin_list$sex + ghsr_exp))
# Significant!
# iii) weight_sac is not linked after accounting for ghsr_exp
anova(lm(weight_sac~ghrelin_list$sex + food_ave + ghsr_exp))
# all significant >:(
# UNSUCCESSFUL

## Try some other stuff!! ## 

## Nhs and Ghrelin expression, shared peak on chr 18
Nhs_geno <- get_genoprob(18, 43.10387)
Nhs_exp <- get_exp_dat("Nhs") #(X)
ghsr_exp <- ghrelin_shortlist$ghrs
triple.fit(Nhs_exp, ghsr_exp, Nhs_geno)
# independent: 2175.568
# reactive: 1757.540
# causal: 1767.019
# complex: 1791.653
# Suggests a reactive model Nhs_geno -> ghsr_exp -> Nhs_exp

# Slight modification so I can use gene expression, food intake, and body weight
#x = gene exp -> pheno1 (foodintake)
# y =pheno -> pheno2 (body weight)
#q = genotype dat -> gene exp (ghsr)
weight_sac <- ghrelin_shortlist$weight_sac
trip <- function(X, Y, Q) {
  bic.independent <- BIC(lm(X~Q)) + BIC(lm(Y~Q)) # X<-Q->Y
  bic.reactive <- BIC(lm(X~Y)) + BIC(lm(Y~Q)) # Q->Y->X
  bic.causal <- BIC(lm(X~Q)) + BIC(lm(Y~X)) # Q->X->Y
  bic.complex <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))
  
  # Print out the scores from each model
  print("BIC Scores of each model")
  scores <- c(bic.independent, bic.reactive, bic.causal, bic.complex)
  names(scores) <- c("independent", "reactive", "causal", "complex")
  print(scores)
}

trip(food_ave, weight_sac, ghsr_exp)
# negative values???!
#independent    reactive      causal     complex 
#-1671.526   -1858.621   -1841.609   -1864.706 
# suggests... complex? 







