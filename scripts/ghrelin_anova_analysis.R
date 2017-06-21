#Week3
#Today I will work more with my ghrelin-related phenotypes and see if I can confirm
#or discover any relationships
#setwd("/Users/s-allens/Documents/ssp/summer_project")

library(tidyverse)
library(plotly)
library(qtl2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)

#loading in mRNA data in order to cut the data set down so you have info from the same animals in each data set
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")

# getting access to scatterplot functions that are stored in a separate file
source("scripts/functions.R")

#load ghrelin phenotype lists
ghrelin_phenos <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_phenos.csv", as.is=TRUE)
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist.csv", as.is=TRUE)

#make use of those ANOVA functions! :)

#test if somatostatin has a significant effect on the phenotypes
#something wrong with Ins_iAUC, will just do shortlist
sst_pvals <- pval_vec_1co(ghrelin_shortlist, "sex", "sst")
sst_sig_lists <- significance(sst_pvals)

#test if food_ave had a significant effect
foodav_pvals <- pval_vec_1co(ghrelin_shortlist, "sex", "food_ave")
foodav_sig_lists <- significance(sst_pvals)

#appetite hormones
#leptin, ghrelin, pyy

#see what affects num islets and insulin per islet
#does number of diet days effect either of these?
#is WPIC (or the other two) related to insulin expression

#get phenotype data
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)

islet_phenos <- matched_phenotypes %>% select(num_islets, Ins_per_islet, WPIC, diet_days)

# find out which phenotypes should be log transformed:
#creates histograms all on one page of the phenotypes 
islet_phenos[,1:4] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")
#these do not look normal, so lets try with a log transform
#creates histograms all on one page of all of the phenotypes (log transformed for normality)
islet_phenos[,1:4] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
  geom_histogram() +
  facet_wrap(~var, scales="free")
#This looks a lot more normal

#convert data into a normal format so more statistical tests can be done
islet_phenos <-log10(islet_phenos[,1:4])
islet_phenos <- cbind(ghrelin_shortlist, islet_phenos)

#correlation matrix
pcor <- cor(islet_phenos[,2:21], use= "complete.obs")
round(pcor, digits=2)
corrplot(pcor, order = "hclust" )
#nothing seems to be strongly correlated with number of islets


#weight growth rate
weights <- matched_phenotypes %>% select(weight_1wk, weight_2wk, weight_3wk, weight_4wk, weight_5wk, weight_6wk, weight_7wk, weight_8wk, weight_9wk, weight_10wk, weight_11wk, weight_12wk, weight_13wk, weight_14wk, weight_15wk, weight_16wk, weight_17wk)
head(weights)

change12 <- weights$weight_2wk - weights$weight_1wk
change23 <- weights$weight_3wk - weights$weight_2wk
change34 <- weights$weight_4wk - weights$weight_3wk
change45 <- weights$weight_5wk - weights$weight_4wk
change56 <- weights$weight_6wk - weights$weight_5wk
change67 <- weights$weight_7wk - weights$weight_6wk
change78 <- weights$weight_8wk - weights$weight_7wk
change89 <- weights$weight_9wk - weights$weight_8wk
change910 <- weights$weight_10wk - weights$weight_9wk
change1011 <- weights$weight_11wk - weights$weight_10wk
change1112 <- weights$weight_12wk - weights$weight_11wk
change1213 <- weights$weight_13wk - weights$weight_12wk
change1314 <- weights$weight_14wk - weights$weight_13wk
change1415 <- weights$weight_15wk - weights$weight_14wk
change1516 <- weights$weight_16wk - weights$weight_15wk
change1617 <- weights$weight_17wk - weights$weight_16wk
weight_change_ave <- (change12 + change23 + change34 + change45 + 
              change56 + change67 + change78 + change89 +
              change910 + change1011 + change1112 + change1213 +
              change1314 + change1415 + change1516 + change1617)/16

islet_phenos <- cbind(islet_phenos, change_ave)




