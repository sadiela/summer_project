# Alpha beta delta phenotype correlations
# Sadie Allen 
# July 21, 2017
# Correlations among cell types and phenotypes
 
rm(list = ls())

library(tidyverse)
library(plotly)
library(qtl2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)

source("scripts/functions.R")

#Load data
#pheno data
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

########################################
# prepare phenotypes
food_6wk <- matched_phenos$food_6wk
food_ave <- matched_phenos$food_ave
weight_6wk <- matched_phenos$weight_6wk
weight_16wk <- matched_phenos$weight_16wk
weight_change <- weight_16wk - weight_6wk

phenos <- data.frame(food_6wk = food_6wk, food_ave = food_ave, weight_6wk = weight_6wk, 
                     weight_16wk = weight_16wk, weight_change = weight_change)

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

#food_6wk, food_ave, weight_16wk, and weight_6wk appear signifcantly more normal after
# log transformation

phenos$food_6wk <- log10(phenos$food_6wk)
phenos$food_ave <- log10(phenos$food_ave)
phenos$weight_6wk <- log10(phenos$weight_6wk)
phenos$weight_16wk <- log10(phenos$weight_16wk)
phenos$weight_change <- phenos$weight_16wk - phenos$weight_6wk

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

# Yay they are all normal!!! :))

########################################

sex <- annot.samples$Sex

# Load in eigengene data
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

gene_frame <- data.frame(alpha = alpha_eigengene, delta = delta_eigengene, beta = Ins2_exp)

gene_pheno <- cbind(gene_frame, phenos)
rownames(gene_pheno) <- annot.samples$Mouse.ID

quartz()
par(mfrow=c(3,2))
ggplot(gene_pheno, aes(x = alpha, y = weight_6wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = alpha, y = weight_16wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = alpha, y = food_ave, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = delta, y = weight_6wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = delta, y = weight_16wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = delta, y = food_ave, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

cor(gene_pheno$alpha, gene_pheno$weight_6wk)
cor(gene_pheno$alpha, gene_pheno$weight_16wk)
cor(gene_pheno$alpha, gene_pheno$food_ave)
cor(gene_pheno$delta, gene_pheno$weight_6wk)
cor(gene_pheno$delta, gene_pheno$weight_16wk)
cor(gene_pheno$delta, gene_pheno$food_ave)

ggplot(gene_pheno, aes(x = beta, y = weight_6wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = beta, y = weight_16wk, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = beta, y = food_ave, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

cor(gene_pheno$beta, gene_pheno$weight_6wk)
cor(gene_pheno$beta, gene_pheno$weight_16wk)
cor(gene_pheno$beta, gene_pheno$food_ave)

ggplot(gene_pheno, aes(x = delta, y = alpha, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = delta, y = beta, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

ggplot(gene_pheno, aes(x = alpha, y = beta, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)

cor(gene_pheno$delta, gene_pheno$alpha)
cor(gene_pheno$delta, gene_pheno$beta)
cor(gene_pheno$alpha, gene_pheno$beta)
































