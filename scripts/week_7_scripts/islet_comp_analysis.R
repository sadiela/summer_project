# Islet Composition Analysis
# Sadie Allen
# July 19, 2017
# Looking at alpha and delta cell eigengenes and genes known to be highly expressed in beta cells

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

food_6wk <- matched_phenos$food_6wk
food_ave <- matched_phenos$food_ave
weight_6wk <- matched_phenos$weight_6wk
weight_16wk <- matched_phenos$weight_16wk
weight_change <- phenos$weight_change

quartz()
phenos[,1:5] %>% 
  gather(var, value) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~var, scales="free")

# Yay they are all normal!!! :))

########################################

# because I do not have an eigengene for beta cells, I will use Ins2 as a surrogate

# Load in eigengene data
alpha_data <- read.csv(file = "data/islet_composition/alpha_eigengene.csv")
alpha_eigengene <- alpha_data$MEgreenyellow

delta_data <- read.csv(file = "data/islet_composition/delta_eigengene.csv")
delta_eigengene <- delta_data$MEyellowgreen

Ins2_exp <- get_exp_dat("Ins2")

gene_frame <- data.frame(alpha = alpha_eigengene, delta = delta_eigengene, beta = Ins2_exp)

gene_pheno <- cbind(gene_frame, phenos)
rownames(gene_pheno) <- annot.samples$Mouse.ID


##################################
### ALPHA EIGENGENE ###
sex <- annot.samples$Sex

# Alpha eigengene peaks: 
QA1 <- get_genoprob(1, 127.21127)
QA6 <- get_genoprob(6, 24.64077)
QA11 <- get_genoprob(11, 16.65394)
QA15 <- get_genoprob(15, 61.93315) 

# Simple BIC Modeling #
BIC(lm(alpha_eigengene~sex))
# -1314.707
BIC(lm(alpha_eigengene~sex + QA1))
# -1302.03
BIC(lm(alpha_eigengene~sex + QA6))
# -1300.296 
BIC(lm(alpha_eigengene~sex + QA11))
# -1301.784
BIC(lm(alpha_eigengene~sex + QA15))
# -1308.781
BIC(lm(alpha_eigengene~sex + weight_16wk))
# -1404.57 
BIC(lm(alpha_eigengene~sex + QA15 + weight_16wk)) 
# -1394.185 

BIC(lm(weight_16wk~sex))
# 2498.531
BIC(lm(weight_16wk~sex + alpha_eigengene))
# 2407.667, significant BIC score drop 
BIC(lm(weight_16wk~sex + QA1))
# 2526.503
BIC(lm(weight_16wk~sex + QA6))
# 2528.206
BIC(lm(weight_16wk~sex + QA11))
# 2518.658
BIC(lm(weight_16wk~sex + QA15))
# 2529.559
# BIC modeling does not appear to work well when you use the eight-way genoprobs. My suspicion
# is that, since you pay a price for the complexity of a model, just using a variable like that
# overrides its effects on the accuracy of the model. I will test this later with the coat color
# example because I know it has a very large peak and the BIC score should definitely drop for it. 

# anova linkage analysis # 
# is the alpha eigengene linked to the body weight phenotypes?
anova(lm(alpha_eigengene~sex + weight_6wk)) # 6.184e-15 ***
anova(lm(alpha_eigengene~sex + weight_16wk)) # < 2.2e-16 ***

# is the alpha eigengene linked to its QTL peaks?
anova(lm(alpha_eigengene~sex + QA1)) # 0.0001786 ***
anova(lm(alpha_eigengene~sex + QA6)) # 0.0003828 ***
anova(lm(alpha_eigengene~sex + QA11)) # 0.0002085 ***
anova(lm(alpha_eigengene~sex + QA15)) # 1.116e-05 ***

# are the body weight phenotypes linked to the alpha eigengene?
anova(lm(weight_6wk~sex + alpha_eigengene)) # 6.184e-15 ***
anova(lm(weight_16wk~sex + alpha_eigengene)) # < 2.2e-16 ***

# are the body weight phenotypes linked to the alpha QTL peaks?
anova(lm(weight_6wk~sex + QA1)) # 0.09196 
anova(lm(weight_6wk~sex + QA6)) # 0.08441
anova(lm(weight_6wk~sex + QA11)) # 0.0004369 *** YES
anova(lm(weight_6wk~sex + QA15)) # 0.1185 

anova(lm(weight_16wk~sex + QA1)) # 0.08947  
anova(lm(weight_16wk~sex + QA6)) # 0.1531
anova(lm(weight_16wk~sex + QA11)) # 0.005431 ** YES
anova(lm(weight_16wk~sex + QA15)) # 0.2284  

# is the link between body weight and the chromosome 11 peak broken after accounting for alpha
# eigengene expression?
anova(lm(weight_6wk~sex +  alpha_eigengene + QA11 )) 
# sex               1 4128.1  4128.1 253.7656 < 2.2e-16 ***
# alpha_eigengene   1 1103.5  1103.5  67.8353 3.142e-15 ***
# QA11              7  269.3    38.5   2.3647    0.0225 *  link not completely broken
anova(lm(weight_16wk~sex +  alpha_eigengene + QA11 )) 
# sex               1  7626.2  7626.2 235.6009 <2e-16 ***
# alpha_eigengene   1  3499.7  3499.7 108.1192 <2e-16 ***
# QA11              7   221.3    31.6   0.9768 0.4478  YES!!

# are weight phenos still linked to alpha eigengene after accounting for QA11?
anova(lm(weight_6wk~sex + QA11 + alpha_eigengene)) # 1.693e-12 *** YES
anova(lm(weight_16wk~sex + QA11 + alpha_eigengene)) # < 2e-16 *** YES

# One successful mediation relationship: QA11 --> alpha_eigengene --> weight_16wk
# Remaining questions: What is under the chromosome 11 peak??

############################################
## DELTA EIGENGENE ##
QD4 <- get_genoprob(4, 13.186538)
QD6 <- get_genoprob(6, 4.776124)
QD13 <- get_genoprob(13, 95.992420)
QD18 <- get_genoprob(18, 5.990333)

# BIC Analysis
BIC(lm(delta_eigengene~sex))
BIC(lm(delta_eigengene~sex + QD6))
# little/no effect
BIC(lm(delta_eigengene~sex + QD18))
# score goes up?!?!

# anova linkage analysis #
# is the delta eigengene linked to the body weight phenotypes?
anova(lm(delta_eigengene~sex + weight_6wk)) # 1.601e-09 ***
anova(lm(delta_eigengene~sex + weight_16wk)) # 1.179e-12 ***

# is the delta eigengene linked to its QTL peaks?
anova(lm(delta_eigengene~sex + QD4)) # 0.0007504 ***
anova(lm(delta_eigengene~sex + QD6)) # 1.223e-07 ***
anova(lm(delta_eigengene~sex + QD13)) # 0.0001694 ***
anova(lm(delta_eigengene~sex + QD18)) # 1.097e-06 ***

# are the body weight phenotypes linked to the delta eigengene?
anova(lm(weight_6wk~sex + delta_eigengene)) # 1.601e-09 ***
anova(lm(weight_16wk~sex + delta_eigengene)) # 1.179e-12 ***

# are the body weight phenotypes linked to the delta QTL peaks?
anova(lm(weight_6wk~sex + QD4)) # 0.2201 
anova(lm(weight_6wk~sex + QD6)) # 0.5467  
anova(lm(weight_6wk~sex + QD13)) # 0.7947 
anova(lm(weight_6wk~sex + QD18)) # 0.9093 

anova(lm(weight_16wk~sex + QD4)) # 0.4549   
anova(lm(weight_16wk~sex + QD6)) # 0.3949 
anova(lm(weight_16wk~sex + QD13)) # 0.5139
anova(lm(weight_16wk~sex + QD18)) # 0.9175 Nopee

# No interesting links between delta eigengene QTLs and body weight

# Hypothesis: body weight of mice affects composition of islets because when mice weigh more,
# they become insulin resistant and are pressured to secrete greater amounts of insulin

########################################
# Alpha and Delta Eigengenes: Peak Locations w Weight Covar
# ALPHA
# lodindex lodcolumn chr      pos      lod    ci_lo     ci_hi
# 1        1     alpha   6 24.64077 6.560362 23.58233 127.26890
# 2        1     alpha  10 28.38945 9.846031 27.51173  28.89694
# 3        1     alpha  15 59.91944 7.582551 58.86843  62.20737
# 4        1     alpha  16 85.93901 6.261213 10.81482  86.91696
# DELTA
#lodindex lodcolumn chr        pos       lod     ci_lo      ci_hi
# 1        1     delta   4  13.186538  6.227576  0.207518  86.964081
# 2        1     delta   6   4.474202 10.261670  0.738935  13.228975
# 3        1     delta  15 101.288524  6.978701 99.774723 102.724056
# 4        1     delta  18   8.012464 10.199439  0.026149   8.530216
##########################################

triple.fit(alpha_eigengene, weight_16wk, QA11)



load("data/qtls_blups_assoc.RData")


rm(list = ls()[grep("^qtl", ls())])













