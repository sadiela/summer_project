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
# Simple BIC Modeling

sex <- annot.samples$Sex

# Alpha eigengene peaks: 
QA1 <- get_genoprob(1, 127.21127)
QA6 <- get_genoprob(6, 24.64077)
QA11 <- get_genoprob(11, 16.65394)
QA15 <- get_genoprob(15, 61.93315) # Why isn't this working???

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

anova(lm(alpha_eigengene~sex))
anova(lm(alpha_eigengene~sex + weight_16wk)) # approx. 0

############################################
anova(lm(alpha_eigengene~sex + QA1)) # 0.0001786 ***  
anova(lm(alpha_eigengene~sex + QA1 + weight_16wk)) 
#               Df  Sum Sq Mean Sq  F value    Pr(>F)    
# sex           1 0.34828 0.34828 261.6248 < 2.2e-16 ***
# QA1           7 0.02844 0.00406   3.0521  8.901e-06 *** 
# weight_16wk   1 0.13338 0.13338 100.1938 < 2.2e-16 ***
# Residuals   368 0.48989 0.00133

anova(lm(alpha_eigengene~sex + QA6)) # 0.0003828 
anova(lm(alpha_eigengene~sex + QA6 + weight_16wk))  
#               Df  Sum Sq Mean Sq F value    Pr(>F)    
# sex           1 0.34828 0.34828 275.224 < 2.2e-16 ***
# QA6           7 0.04514 0.00645   5.096 1.529e-05 ***
# weight_16wk   1 0.14089 0.14089 111.335 < 2.2e-16 ***
# Residuals   368 0.46569 0.00127                      

anova(lm(alpha_eigengene~sex + QA11)) # 0.0002085 ***
anova(lm(alpha_eigengene~sex + QA15)) # 1.116e-05 ***
############################################

anova(lm(weight_16wk~sex))

############## Successful mediation....?
anova(lm(weight_16wk~sex + alpha_eigengene)) # approx. 0
anova(lm(weight_16wk~sex + QA1)) # 0.0322, weak, but there....
#             Df  Sum Sq Mean Sq  F value  Pr(>F)    
# sex         1  7626.2  7626.2 187.5826 < 2e-16 ***
# QA1         7   631.1    90.2   2.2177 0.08947 .  
# Residuals 369 15001.8    40.7  
anova(lm(weight_16wk~sex + QA1 + alpha_eigengene)) # All significant!?!
#                   Df  Sum Sq Mean Sq  F value    Pr(>F)    
# sex               1  7626.2  7626.2 238.0081 < 2.2e-16 ***
# QA1               7   631.1    90.2   2.8139  0.03024 * 
# alpha_eigengene   1  3210.4  3210.4 100.1938 < 2.2e-16 ***
# Residuals       368 11791.4    32.0  
anova(lm(weight_16wk~sex + alpha_eigengene + QA1)) # All significant!?!
#########################################


anova(lm(weight_16wk~sex + QA6)) # 0.1531


##### SUCCESSFUL MEDIATION ####
anova(lm(weight_16wk~sex + QA11)) # 0.005431
anova(lm(weight_16wk~sex + QA11 + alpha_eigengene)) # 0.005431
#                   Df  Sum Sq Mean Sq  F value  Pr(>F)    
# sex               1  7626.2  7626.2 235.6009 < 2e-16 ***
# QA11              7   822.0   117.4   3.6277 0.00085 ***
# alpha_eigengene   1  2899.1  2899.1  89.5628 < 2e-16 ***
# Residuals       368 11911.9    32.4 
anova(lm(weight_16wk~sex +  alpha_eigengene + QA11 )) # 0.005431
# OMG A MEDIATOR!!!! :))))) YAYYYYYAYAYAY
anova(lm(weight_16wk~sex + QA15)) # 0.2284


# Delta eigengene peaks:
QD4 <- get_genoprob(4, 13.186538)
QD6 <- get_genoprob(6, 4.776124)
QD13 <- get_genoprob(13, 95.992420)
QD18 <- get_genoprob(18, 5.990333)

BIC(lm(delta_eigengene~sex))
BIC(lm(delta_eigengene~sex + QD6))
# little/no effect
BIC(lm(delta_eigengene~sex + QD18))
# score goes up?!?!

anova(lm(delta_eigengene~sex))
anova(lm(delta_eigengene~sex + QD6)) # significant
anova(lm(delta_eigengene~sex + weight_16wk)) # significant
anova(lm(delta_eigengene~sex + weight_16wk + QD4)) # all significant
anova(lm(delta_eigengene~sex + weight_16wk + QD6)) # all significant
anova(lm(delta_eigengene~sex + weight_16wk + QD13)) # all significant
anova(lm(delta_eigengene~sex + weight_16wk + QD18)) # all significant

anova(lm(weight_16wk~sex))
anova(lm(weight_16wk~sex + delta_eigengene)) #signif
anova(lm(weight_16wk~sex + alpha_eigengene)) #signif
anova(lm(weight_16wk~sex + Ins2_exp)) #signif
anova(lm(weight_16wk~sex + delta_eigengene + alpha_eigengene + Ins2_exp)) #all significant

anova(lm(weight_16wk~sex + QD4)) # 0.4549
anova(lm(weight_16wk~sex + QD6)) # 0.3949
anova(lm(weight_16wk~sex + QD13)) # 0.5139
anova(lm(weight_16wk~sex + QD18)) # 0.9175



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




















