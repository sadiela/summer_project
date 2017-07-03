# Exploring Genes of Interest: Ghsr, Food_ave, Weight_sac
# Sadie Allen
# June 29, 2017
# In this script I will continue my investigation into the genes of interest identified
# in my p-value tables and combine this information with the location of the peak on chr 18
# for Ghsr expression

rm(list = ls())

# Load libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)

# Load all functions
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


# find genes under the peak, compare with significant list of genes for Ghsr
# FINDING WHAT'S UNDER THE PEAK

# Load in ensembl genes
peak_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/chr_18_ghsr_gup.csv")

# Choose all protein-coding genes
peak_genes_pc <- peak_genes[peak_genes$Gene.type == "protein_coding", ]

# save to a new data file
write.csv(peak_genes_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/chr_18_ghsr_gup_pc.csv")
peak_genes_pc <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/chr_18_ghsr_gup_pc.csv")

ghsr_gene_peak_names <- as.character(peak_genes_pc$Gene.name)

# Get list of genes significant to Ghsr expression
# Calculate pvalues for all gene expressions
ghsr_pvals <- pvals_changing_covar(ghrelin_list$ghsr_exp, rankz.mrna, ghrelin_list$sex)

# Get gene names
gene_names <- get_gene_names()

# set names
names(ghsr_pvals) <- gene_names

# Get rid of insignificant genes and sort by p-value
ghsr_pvals <- sig_list(ghsr_pvals, .0001)
head(ghsr_pvals, n = 10)
#write.csv(ghsr_pvals, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghsr_pvals.csv")
#ghsr_pvals <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghsr_pvals.csv")

sig_pval_ghsr <- names(ghsr_pvals)
head(sig_pval_ghsr)
head(ghsr_gene_peak_names)

candidate_genes <- intersect(sig_pval_ghsr, ghsr_gene_peak_names)

# Do qtl scans of candidate genes to see which ones have local peaks on chromosome 18
# Create data frame of expression data for all of the candidate genes
cand_exp_data <- data.frame(mouse.id = ghrelin_list$Mouse.ID)
for(i in 1:length(candidate_genes)) {
  expr_data <- get_exp_dat(candidate_genes[i])
  cand_exp_data <- cbind(cand_exp_data, expr_data)
}
names <- c("mouse.id", candidate_genes)
colnames(cand_exp_data) <- names

# QTL scan prep

load("data/qtl_prep.RData")

#probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Calculate kinship
#kin <- calc_kinship(probs = probs, type = "loco", cores = 4)

# Additive covariates
#temp = merge(annot.samples, ghrelin_list, by = "row.names")
#rownames(temp) = temp[,1]
# annot.samples <- merge(annot.samples, diet_days_id, by.x="Mouse.ID", by.y = "Mouse.ID")
# Now, create the covariate
#add_covar <- model.matrix(~Sex + Generation + diet_days, data = temp)[,-1]

# Convert a marker map organized as data frame to a list
#map <- map_df_to_list(map = snps, pos_column = "bp")

# Save all this stuff so you can just load it in!!
#save(probs, snps, kin, add_covar, map, file = "data/qtl_prep.RData")
#load("data/qtl_prep.RData")
qtl.Ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[, colnames(ghrelin_list) == "ghsr_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Thoc1 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Thoc1", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Ccny <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Ccny", drop = FALSE],
               kinship = kin, addcovar = add_covar, cores = 4)
qtl.Cul2 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Cul2", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Mtpap <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Mtpap", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Usp14 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Usp14", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Rpl27ps3 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Rpl27-ps3", drop = FALSE],
                      kinship = kin, addcovar = add_covar, cores = 4)
qtl.Fzd8 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Fzd8", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Svil <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Svil", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Mpp7 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Mpp7", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Epc1 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Epc1", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
qtl.Rock1 <- scan1(genoprobs = probs, pheno = cand_exp_data[,colnames(cand_exp_data)=="Rock1", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)

qtl_df <- data.frame(Thoc1 = qtl.Thoc1, Ccny = qtl.Ccny, Cul2 = qtl.Cul2, Mtpap = qtl.Mtpap, 
                     Usp14 = qtl.Usp14, Rpl27ps3 = qtl.Rpl27ps3, Fzd8 = qtl.Fzd8, Svil = qtl.Svil, 
                     Mpp7 = qtl.Mpp7, Epc1 = qtl.Epc1, Rock1 = qtl.Rock1) 

qtl_18df <- qtl_df[grep("^18_", rownames(qtl_df)),]
max_lod <- apply(qtl_18df, 2, max)
max_lod_position <- apply(qtl_18df, 2, which.max)

remaining_candidates <- names(max_lod[max_lod > 6])

# Plot remaining candidates
quartz()
par(mfrow = c(3, 1))
plot_scan1(x = qtl.Ccny, map = map, main = colnames(qtl.Ccny)[1])
plot_scan1(x = qtl.Mtpap, map = map, main = colnames(qtl.Mtpap)[1])
plot_scan1(x = qtl.Usp14, map = map, main = colnames(qtl.Usp14)[1])
quartz()
par(mfrow = c(3, 1))
plot_scan1(x = qtl.Fzd8, map = map, main = colnames(qtl.Fzd8)[1])
plot_scan1(x = qtl.Svil, map = map, main = colnames(qtl.Svil)[1])
plot_scan1(x = qtl.Mpp7, map = map, main = colnames(qtl.Mpp7)[1])

quartz()
par(mfrow =c(2, 1))
plot_scan1(x = qtl.Ghsr, map = map, chr = 18, main = colnames(qtl.Ghsr)[1])
plot_scan1(x = qtl.Svil, map = map, chr = 18, main = colnames(qtl.Svil)[1])


find_peaks(qtl.Ccny, map, threshold = 6, drop = 1.5) 
# pos: 9.98, lod = 10.24
find_peaks(qtl.Mtpap, map, threshold = 6, drop = 1.5) 
# pos: 4.303, lod = 16.049
find_peaks(qtl.Usp14, map, threshold = 6, drop = 1.5) 
# pos 83.95, lod = 6 ~~~~~ THIS ONE IS IFFY
find_peaks(qtl.Fzd8, map, threshold = 6, drop = 1.5) 
# pos 6.6, lod = 10.4
find_peaks(qtl.Svil, map, threshold = 6, drop = 1.5) 
# pos 4.77, lod = 49.236
find_peaks(qtl.Mpp7, map, threshold = 6, drop = 1.5) 
# pos 7.00, lod = 32.8

# GR Interval:  (interval 0.26149 - 10.21129)
# Therefore, I can eliminate Usp14 as a candidate gene
remaining_candidates <- remaining_candidates[-3]

# Get genotype as the positions of the peaks for all of the remaining candidates
Ccny_geno <- get_genoprob(18, 9.983426)
Mtpap_geno <- get_genoprob(18, 4.303864)
Fzd8_geno <- get_genoprob(18, 6.657569)
Svil_geno <- get_genoprob(18, 4.772599)
Mpp7_geno <- get_genoprob(18, 7.003513)


# BIC Modeling! :) 
# Remember, triple.fit takes: 
# X: gene expression data for a given gene
# Y: quantitative measurement of clinical phenotype
# Q: genotype at a given marker

# ghsr expression will function as the clin phenotype(Y), gene expression of the gene as X, and genotype as Q
triple.fit(cand_exp_data$Ccny, ghrelin_list$ghsr_exp, Ccny_geno)
# Suggests a reactive model: Q -> Y -> X (genotype -> ghsr_exp -> ccny_exp)

triple.fit(cand_exp_data$Mtpap, ghrelin_list$ghsr_exp, Mtpap_geno)
# inconclusive, (reactive - causal < 5)

triple.fit(cand_exp_data$Fzd8, ghrelin_list$ghsr_exp, Fzd8_geno)
# suggests a causal model: Q -> X -> Y (genotype -> Fzd8_exp -> ghsr_exp)

triple.fit(cand_exp_data$Svil, ghrelin_list$ghsr_exp, Svil_geno)
# suggests a causal model: Q -> X -> Y (genotype -> Svil_exp -->  ghsr_exp)

triple.fit(cand_exp_data$Mpp7, ghrelin_list$ghsr_exp, Mpp7_geno)
# inconclusive, (reactive - causal < 5)


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

###########
# Trying to confirm svil_expr mediates the effect of the chr18 peak genotype on ghsr expression
Q18 <- get_genoprob(18, 4.772599)
svil_exp <- get_exp_dat("Svil")
ghsr_exp <- ghrelin_list$ghsr_exp
sex <- ghrelin_list$sex
# i) ghsr expression is linked to Q18
anova(lm(ghsr_exp~sex + Q18))
# significant (***)
# ii) svil_exp is linked to Q18
anova(lm(svil_exp~sex + Q18))
# significant (***)
# iii) ghsr_exp not linked after accounting for Q18
anova(lm(ghsr_exp~sex + svil_exp + Q18))
# kindd of significant? (*)
# iv) svil gene expression still linked after accounting for ghsr expression
anova(lm(svil_exp~sex + ghsr_exp + Q18))
# significant (***)

# Trying to confirm Mpp7_exp mediates the effect of the chr18 peak genotype on ghsr expression
Q18 <- get_genoprob(18, 7.003513)
Mpp7_exp <- get_exp_dat("Mpp7")
# i) ghsr expression is linked to Q18
anova(lm(ghsr_exp~sex + Q18))
# Genotype at this locus not very significantly linked to ghsr expression... use a slightly 
# different locus? 

# Trying to confirm that ghsr_expr mediates the effect of the genotype (Ccny_geno) of the ccny peak on ccny expression Ccny_exp
# i) Ccny expression is linked to Ccny_geno
anova(lm(cand_exp_data$Ccny~ghrelin_list$sex + Ccny_geno))
# not significant!?!?!?

# Trying to confirm that Fzd8 expression mediates the effect of the Fzd8 peak genotype on ghsr expression
# i) Ghsr expression is linked to Fzd8_geno
anova(lm(ghrelin_list$ghsr_exp~ghrelin_list$sex + Fzd8_geno))
# significant (**)
# ii) Fzd8 expression is linked to Fzd8_geno
anova(lm(cand_exp_data$Fzd8~ghrelin_list$sex + Fzd8_geno))
# significant (***)
# iii) Ghsr_exp not linked after accounting for Fzd8_exp
anova(lm(ghrelin_list$ghsr_exp~ghrelin_list$sex + cand_exp_data$Fzd8+ Fzd8_geno))
# All are still significant >:(







