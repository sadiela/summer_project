# GOI/POI peak genes analysis
# Sadie Allen
# July 11, 2017
# Analyzing and narrowing down genes under the peaks of significant genes and phenotypes

rm(list = ls())

library(corrplot)

ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

source("scripts/functions.R")

## GENES UNDER PEAKS ##
# STILL HAVEN'T DONE ANYTHING WITH THESE
# Efnb3 chromosome 4:
# Trp53inp1, Rad54b, Triqk, Plekhf2

# Efnb3 chromosome 18:
# Map3k8

# Ghsr chromosome 4:
# Rbm12b2, 1110037F02Rik, Plag1, Ints8, Rbm12b1, Ccne2, Tgs1, Clvs1, Nsmaf, Gm11808, Ubxn2b, Pdp1, Cdh17, Impad1, Asph

# Ghsr chromosome 18:
# Thoc1, Ccny, Cul2, Mtpap, Usp14, Rpl27-ps3, Fzd8, Svil

# food_ave chromosome 7: 
# Mki67



### GENE PEAK MEDIATORS ###
# Mediators of Efnb3 chr 4 peak:
# Ghsr, Sst, Cacna1h
# Mediators of Efnb3 chr 18 peak:
# Nhs, Cacna1h, Crhr2, Ghsr
# Mediators of Ghsr Chr 18 peak:
# Nhs, Cacna1h, Kcnd2, Efnb3
# Mediators of Ghsr chr 4 peak:
# Arg1

ghsr_exp <- get_exp_dat("Ghsr")
efnb3_exp <- get_exp_dat("Efnb3")
sst_exp <- get_exp_dat("Sst")
cacna1h_exp <- get_exp_dat("Cacna1h")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
kcnd2_exp <- get_exp_dat("Sst")
arg1_exp <- get_exp_dat("Arg1")
food_ave <- ghrelin_list$food_ave
weight_sac <- ghrelin_list$weight_sac
weight_change <- ghrelin_list$weight_change
weight_6wk <- ghrelin_list$weight_6wk

mediator_genes <- data.frame(row.names = rownames(ghrelin_list), ghsr = ghsr_exp, efnb3 = efnb3_exp,
                             sst = sst_exp, cacna1h = cacna1h_exp, nhs = nhs_exp, crhr2 = crhr2_exp,
                             kcnd2 = kcnd2_exp, arg1 = arg1_exp, food_ave = food_ave, weight_sac = weight_sac,
                             weight_change = weight_change, weight_6wk = weight_6wk)

# Make a correlation plot! :P
pcor <- cor(mediator_genes, use= "complete.obs")
pcor <- round(pcor, digits = 2)
corrplot(pcor, order = "hclust")


# Arg1 does not have Efnb3 or Ghsr as mediators
# Sst has a peak on 6 and 18


# Plug these into mediation analysis to verify it 

# reverse covariate scans to see who's influencing who
# the gene that drops the other gene's peak is the upstream gene



