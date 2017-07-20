# Week2
# scatterplots and correlation matrices of ghrelin pathway-related phenotypes
# setwd("/Users/s-allens/Documents/ssp/summer_project")

rm(list = ls())

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


#This section works
#saving matched phenotypes to a new data file
#write.csv(matched_phenotypes, file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", row.names = FALSE)
#NOW THAT THIS HAS BEEN DONE ONCE, SIMPLY RUN THE FOLLOWING TO GET ACCESS TO THE MATCHED PHENOTYPES
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/old_data/matched_pheno_clin.csv", as.is=TRUE)

scatter_phenotypes <- matched_phenotypes %>% select(num_islets, Ins_per_islet, WPIC, Glu_0min, 
                      Glu_6wk, Glu_10wk, Glu_14wk, Ins_0min, Ins_6wk, Ins_10wk,
                      Ins_14wk, Glu_tAUC, Glu_iAUC, Ins_tAUC, Ins_iAUC, TG_6wk,
                      TG_10wk, TG_14wk, oGTT_weight, Glu_sac, Ins_sac, TG_sac)

condensed_scatter_phenos <- matched_phenotypes %>% select(num_islets, Ins_per_islet,
                            WPIC, Glu_0min, Glu_14wk, Ins_0min, Ins_14wk, Glu_tAUC,
                            Glu_iAUC, Ins_tAUC, Ins_iAUC, oGTT_weight)

# plot phenotypes
quartz()
pairs(scatter_phenotypes[1:22], upper.panel=panel.cor,diag.panel=panel.hist)

#I looked at all of the correlation plots and deemed which ones looked most interesting.
#This consisted mostly of removing intuitive data, such as the correlation between Ins_6wk
#and Ins_10wk. I also eliminated phenotypes that did not have any strong correlations with other
#phenotypes
quartz()
pairs(condensed_scatter_phenos, upper.panel=panel.cor,diag.panel=panel.hist)

#phenotypes that will be relevant to my ghrelin study
ghrelin_phenos <- matched_phenotypes %>% select(sex, Glu_0min, Glu_tAUC, Glu_iAUC, Glu_6wk, 
                  Glu_10wk, Glu_14wk, Glu_sac, Ins_0min, Ins_tAUC, Ins_iAUC, Ins_6wk, 
                  Ins_10wk, Ins_14wk, Ins_sac, food_ave, weight_sac, G33_ins_secrete, 
                  G83_ins_secrete, G167_ins_secrete, KCl_G33_ins_secrete, 
                  GLP1_G83_ins_secrete, AA_G83_ins_secrete, PA_G167_ins_secrete)
ghrelin_shortlist <- matched_phenotypes %>% select(sex, Glu_0min,  Glu_sac, Ins_0min, 
                    Ins_sac, food_ave, weight_sac, G33_ins_secrete, 
                    G83_ins_secrete, G167_ins_secrete)

# find out which phenotypes should be log transformed:
#creates histograms all on one page of all of the phenotypes 
ghrelin_phenos[,2:24] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")
#these do not look normal, so lets try with a log transform
#creates histograms all on one page of all of the phenotypes (log transformed for normality)
ghrelin_phenos[,2:24] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
  geom_histogram() +
  facet_wrap(~var, scales="free")
#This looks a lot more normal

ghrelin_shortlist[,6:8] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

ghrelin_shortlist[,6:8] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
  geom_histogram() +
  facet_wrap(~var, scales="free")


#convert data into a normal format so more statistical tests can be done
sex <- ghrelin_phenos[,1]
ghrelin_phenos <- cbind(sex, log10(ghrelin_phenos[,2:24]))
ghrelin_shortlist <- cbind(sex, log10(ghrelin_shortlist[,2:10]))

pairs(ghrelin_shortlist[2:10], upper.panel=panel.cor, diag.panel=panel.hist)
#be cautious using the Ins_0min phenotype, data is a little funky

#Now, how do I use the mRNA data?

#retrieve expression data for all genes of interest
ghrl_exp <- get_exp_dat("Ghrl")
ghsr_exp <- get_exp_dat("Ghsr")
sst_exp <- get_exp_dat("Sst")
ins1_exp <- get_exp_dat("Ins1")
ins2_exp <- get_exp_dat("Ins2")
sstr3_exp <- get_exp_dat("Sstr3")

#create a dataframe with expression data
gene_expressions <- data.frame(ghrl = ghrl_exp, ghsr = ghsr_exp, 
                     sst = sst_exp, ins1 = ins1_exp, 
                     ins2 = ins2_exp, sstr = sstr3_exp)

#combine with phenotype data
ghrelin_phenos <- cbind(ghrelin_phenos, gene_expressions)
ghrelin_shortlist <- cbind(ghrelin_shortlist, gene_expressions)

#save these data frames since you use them a lot
write.csv(ghrelin_phenos, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_phenos.csv", row.names = FALSE)
write.csv(ghrelin_shortlist, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist.csv", row.names = FALSE)
#NOW THAT THIS HAS BEEN DONE ONCE, SIMPLY RUN THE FOLLOWING TO GET ACCESS TO THE MATCHED PHENOTYPES
ghrelin_phenos <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_phenos2.csv", as.is=TRUE)
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv", as.is=TRUE)

pairs(gene_expressions[1:6], upper.panel = panel.cor, diag.panel = panel.hist)
pairs(ghrelin_shortlist, upper.panel = panel.cor, diag.panel = panel.hist)


#correlation matrix
pcor <- cor(ghrelin_shortlist[2:], use= "complete.obs")
round(pcor, digits=2)
corrplot(pcor, order = "hclust" )


Mouse.ID <- matched_phenotypes$Mouse.ID

ghrelin_phenos <- cbind(Mouse.ID, ghrelin_phenos)
ghrelin_shortlist <- cbind(Mouse.ID, ghrelin_shortlist)

write.csv(ghrelin_phenos, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_phenos2.csv", row.names = FALSE)
write.csv(ghrelin_shortlist, file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv", row.names = FALSE)

