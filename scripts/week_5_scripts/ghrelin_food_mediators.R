# Mediators of Ghrelin Receptor - Food Intake Relationship
# Sadie Allen
# July 2, 2017
# In this script I will look for genes that may mediate the relationship between 
# Ghsr expression and food intake

# clear environment
rm(list = ls())

# load libraries
library(ggplot2)
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)
library(dplyr)

# load functions
source("scripts/functions.R")

# Load data
# islet data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# UPDATED PHENOTYPE DATA
matched_phenos <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv", as.is=TRUE)
rownames(matched_phenos) <- matched_phenos[,1]
ghrelin_list <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list[,2]

# Get list of significant genes for ghsr and food_ave
ghsr_pvals <- pvals_changing_covar(ghrelin_list$ghsr_exp, rankz.mrna, ghrelin_list$sex)

pheno_pvals <- read.csv("results/prelim_research/pheno_gene_exp_signif.csv")
food_pvals <- pheno_pvals$food_ave

gene_names <- get_gene_names()

names(ghsr_pvals) <- gene_names
names(food_pvals) <- gene_names

# Narrow these down so they are just the significant genes
ghsr_pvals <- sig_list(ghsr_pvals, .00001)
food_pvals <- sig_list(food_pvals, .00001)

sig_ghsr_genes <- names(ghsr_pvals)
sig_food_pvals <- names(food_pvals)

# Find genes related to both 
candidate_genes <- intersect(sig_food_pvals, sig_ghsr_genes)
# 586 candidates... well

# get rid of genes that aren't highly correlated with both ghsr_exp and food_ave
narrowed_cand <- character()
for(i in 1:length(candidate_genes)) {
  gene_name <- candidate_genes[i]
  gene_exp <- get_exp_dat(gene_name)
  if(abs(cor(gene_exp, ghrelin_list$ghsr_exp)) > 0.35 && abs(cor(gene_exp, ghrelin_list$food_ave)) > 0.35) {
    narrowed_cand <- c(narrowed_cand, gene_name)
  }
}


sex <- ghrelin_list$sex
food_ave <- ghrelin_list$food_ave
ghsr_exp <- ghrelin_list$ghsr_exp
# Create a mediation analysis loop!
# First: mediation analysis function
# Tests if middle mediates the relationship between begin and end
mediation_analysis <- function(begin, middle_name, end) {
  middle <- get_exp_dat(middle_name)
  # begin is linked to end
  step1_val <- anova(lm(end ~ sex + begin))[2,5]
  # middle is linked to begin
  step2_val <- anova(lm(middle ~ sex + begin))[2,5]
  # end NOT linked to begin after accounting for middle
  step3_val <- anova(lm(end ~ sex + middle + begin))[3,5]
  # middle and begin linked after accounting for end
  step4_val <- anova(lm(middle ~ sex + end + begin))[3,5]
  if(step1_val < 0.001 && step2_val < 0.001 && step3_val > 0.01 && step4_val < 0.001) {
    print(paste(middle_name, "is a mediator of the relationship between ghsr expression and food intake", sep = " "))
    return(TRUE)
  } else {
    print(paste(middle_name, "is NOT a mediator", sep = " "))
    return(FALSE)
  }
}

# Now create a loop
mediators <- character()
for(i in 1:length(narrowed_cand)) {
  print(narrowed_cand[i])
  if(mediation_analysis(ghsr_exp, narrowed_cand[i], food_ave) == TRUE) {
    mediators <- c(mediators, narrowed_cand[i])
  }
}

# This gives us one mediator: Efnb3
# Retest manually

Efnb3 <- get_exp_dat("Efnb3")
anova(lm(food_ave ~ sex + ghsr_exp))
# middle is linked to begin
anova(lm(Efnb3 ~ sex + ghsr_exp))
# end NOT linked to begin after accounting for middle
anova(lm(food_ave ~ sex + Efnb3 + ghsr_exp)) # a weaker link
# middle and begin linked after accounting for end
anova(lm(Efnb3 ~ sex + food_ave + ghsr_exp))



