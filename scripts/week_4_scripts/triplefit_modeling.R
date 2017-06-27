# Attempting linear models
# Sadie Allen
# June 27, 2017
# In this script I will be trying to adapt the triplefit function used on the BTBR 
# f2 cross for use in the DO dataset

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
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load all functions
source("scripts/functions.R")

## Ghrelin & Food Intake ##
genotype_ghsr <- get_genoprob("18_4474313")
ghsr_exp <- ghrelin_shortlist$ghsr
food_ave <- ghrelin_shortlist$food_ave

triple.fit(ghsr_exp, food_ave, genotype_ghsr)
# independent: 175.06208
# reactive: 71.45350
# causal: 47.18930 Suggests a causal relationship! genotype -> gene expression -> phenotype
# complex: 77.67778

## Try some other stuff!! ## 

## Nhs and Ghrelin expression, shared peak on chr 18
Nhs_geno <- get_genoprob(18, 43.10387)
Nhs_exp <- get_exp_dat("Nhs") #(X)
ghsr_exp <- ghrelin_shortlist$ghrs
triple.fit(Nhs_exp, ghsr_exp, Nhs_geno)
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








