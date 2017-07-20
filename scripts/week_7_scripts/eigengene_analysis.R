# Eigengene Analysis
# Sadie Allen
# July 19, 2017
# Investigation of the yellowgreen module created by Petr and Mark that seems to be representative
# of delta cells in the islet tissue of the mice

rm(list = ls())

# Load data
# Yellowgreen module
yellowgreen <- read.csv(file = "data/yellowgreen_rnaseq_module.csv")
# Phenotype Data
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

diet_days <- matched_phenos$diet_days

# Load libraries
library(ggplot2)
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)
library(dbplyr)

# Get list of numeric phenotypes
numeric_phenos <- matched_phenos
x <- vector("numeric")
for(i in 1:length(numeric_phenos)){
  if(is.numeric(numeric_phenos[,i])== TRUE && all(is.na(numeric_phenos[,i])) == FALSE){
    x <- c(x,i)
  }
}
numeric_phenos <- numeric_phenos[,x]
numeric_phenos <- numeric_phenos[, -grep("DOwave", colnames(numeric_phenos))]
numeric_phenos <- numeric_phenos[,-1]


# Test correlations with ALL phenotypes
x <- numeric()
corvals <- numeric()
for(i in 1:ncol(numeric_phenos)) {
  if(all(is.na(numeric_phenos[,i])) == FALSE) {
    if(abs(cor(yellowgreen$MEyellowgreen, numeric_phenos[,i], use = "complete.obs")) > 0.30) {
      corr <- cor(yellowgreen$MEyellowgreen, numeric_phenos[,i], use = "complete.obs")
      corvals <- c(corvals, corr)
      x <- c(x, i)
    }
  }
}
correlated_phenos <- numeric_phenos[,x]
colnames(correlated_phenos)
names(corvals) <- colnames(correlated_phenos)
sort(corvals, decreasing = TRUE)
# a positive correlation with glucagon content, strongest (negative) correlations with weight phenotypes

# Correlated with glucose, insulin, food, and weight-related phenotypes
# Strongest correlations with weight-related phenotypes

#QTLs for weight_20wk, food_ave, weight_change, weight_6wk

module_phenos <- correlated_phenos[,c(18, 33, 44, 54)]

ygdata <- yellowgreen$MEyellowgreen
weight_change <- correlated_phenos$weight_16wk - correlated_phenos$weight_1wk
yellowgreen_phenos <- cbind(module_phenos, weight_change, ygdata)
yellowgreen_phenos$Mouse.ID <- NULL

pcor <- cor(yellowgreen_phenos, use= "complete.obs")
round(pcor, digits=2)
corrplot(pcor, order = "hclust")

write.csv(yellowgreen_phenos, file = "data/yellowgreen_phenos.csv")
yellowgreen_phenos <- read.csv(file = "data/yellowgreen_phenos.csv")


# QTL scan prep
load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

qtl.ygdata <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "ygdata", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.ygdata, map = map, threshold = 6, drop = 1.5)

qtl.food_6wk <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "food_6wk", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_6wk, map = map, threshold = 6, drop = 1.5)

qtl.food_ave <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "food_ave", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.food_ave, map = map, threshold = 6, drop = 1.5)

qtl.weight_6wk <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_6wk", drop = FALSE],
                    kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_6wk, map = map, threshold = 6, drop = 1.5)

qtl.weight_16wk <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_16wk", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_16wk, map = map, threshold = 6, drop = 1.5)

qtl.weight_change <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_change", drop = FALSE],
                        kinship = kin, addcovar = add_covar, cores = 4)
find_peaks(qtl.weight_change, map = map, threshold = 6, drop = 1.5)

quartz()
par(mfrow = c(3,2))
plot_scan1(x = qtl.ygdata, map = map, main = "Yellow Green Module", col = "black")
plot_scan1(x = qtl.food_6wk, map = map, main = "Food 6 Weeks", col = "black")
plot_scan1(x = qtl.food_ave, map = map, main = "Food Average", col = "black")
plot_scan1(x = qtl.weight_6wk, map = map, main = "Weight 6 Weeks", col = "black")
plot_scan1(x = qtl.weight_16wk, map = map, main = "Weight 16 Weeks", col = "black")
plot_scan1(x = qtl.weight_change, map = map, main = "Weight Change", col = "black")

# covariate scans for the phenotypes
temp = merge(annot.samples, yellowgreen_phenos, by = "row.names")
rownames(temp) = temp[,1]
# Now, create the covariates

add_covar_ygdata <- model.matrix(~Sex + Generation + diet_days + ygdata, data = temp)[,-1]

qtl.food6wk_ygdatacovar <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "food_6wk", drop = FALSE],
                             kinship = kin, addcovar = add_covar_ygdata, cores = 4)
qtl.foodave_ygdatacovar <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "food_ave", drop = FALSE],
                                 kinship = kin, addcovar = add_covar_ygdata, cores = 4)
qtl.weight6wk_ygdatacovar <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_6wk", drop = FALSE],
                                   kinship = kin, addcovar = add_covar_ygdata, cores = 4)
qtl.weight16wk_ygdatacovar <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_16wk", drop = FALSE],
                                   kinship = kin, addcovar = add_covar_ygdata, cores = 4)
qtl.weightchange_ygdatacovar <- scan1(genoprobs = probs, pheno = yellowgreen_phenos[, colnames(yellowgreen_phenos) == "weight_change", drop = FALSE],
                                    kinship = kin, addcovar = add_covar_ygdata, cores = 4)

quartz()
plot_scan1(x = qtl.food_6wk, map = map, main = "food6wk vs ygdata covar", col = "black")
plot_scan1(x = qtl.food6wk_ygdatacovar, map = map, add = TRUE, col = "red")

quartz()
plot_scan1(x = qtl.food_ave, map = map, main = "foodave vs ygdata covar", col = "black")
plot_scan1(x = qtl.foodave_ygdatacovar, map = map, add = TRUE, col = "red")

quartz()
plot_scan1(x = qtl.weight_6wk, map = map, main = "weight6wk vs ygdata covar", col = "black")
plot_scan1(x = qtl.weight6wk_ygdatacovar, map = map, add = TRUE, col = "red")

quartz()
plot_scan1(x = qtl.weight_16wk, map = map, main = "weight16wk vs ygdata covar", col = "black")
plot_scan1(x = qtl.weight16wk_ygdatacovar, map = map, add = TRUE, col = "red")

quartz()
plot_scan1(x = qtl.weight_change, map = map, main = "weightchange vs ygdata covar", col = "black")
plot_scan1(x = qtl.weightchange_ygdatacovar, map = map, add = TRUE, col = "red")

# No noticeable effects on peaks



