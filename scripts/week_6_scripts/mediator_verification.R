# Mediator Verification
# Sadie Allen
# July 11, 2017
# Mediation analysis to verify the mediators identified by looking for QTL drops

rm(list = ls())

#Load Libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
library(RSQLite)
library(dbplyr)

# Load data
# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load functions
source("scripts/functions.R")

sex <- annot.samples$Sex

mediation_analysis <- function(begin, middle, end) {
  # begin is linked to end
  step1_val <- anova(lm(end ~ sex + begin))[2,5]
  # middle is linked to begin
  step2_val <- anova(lm(middle ~ sex + begin))[2,5]
  # end NOT linked to begin after accounting for middle
  step3_val <- anova(lm(end ~ sex + middle + begin))[3,5]
  # middle and begin linked after accounting for end
  step4_val <- anova(lm(middle ~ sex + end + begin))[3,5]
  if(step1_val < 0.01 && step2_val < 0.01 && step3_val > 0.05 && step4_val < 0.01) {
    print("is a mediator!!!!!!!!!!!!!")
    return(TRUE)
  } else {
    print("is NOT a mediator")
    return(FALSE)
  }
}

# Get genoprobs at 4 peaks
efnb3_4 <- get_genoprob(4, 5.66854)
efnb3_18 <- get_genoprob(18, 3.351)
ghsr_4 <- get_genoprob(4, 11.984686)
ghsr_18 <- get_genoprob(18, 4.474313)


# Get gene expressions
ghsr_exp <- get_exp_dat("Ghsr")
efnb3_exp <- get_exp_dat("Efnb3")
sst_exp <- get_exp_dat("Sst")
cacna1h_exp <- get_exp_dat("Cacna1h")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
kcnd2_exp <- get_exp_dat("Sst")
arg1_exp <- get_exp_dat("Arg1")
svil_exp <- get_exp_dat("Svil")


# Relationships to test:
#efnb3_4 -> Ghsr -> efnb3_exp
mediation_analysis(efnb3_4, "Ghsr", efnb3_exp)
#efnb3_4 -> sst -> efnb3_exp
mediation_analysis(efnb3_4, "Sst", efnb3_exp)
#efnb3_4 -> cacna1h -> efnb3_exp
mediation_analysis(efnb3_4, "Cacna1h", efnb3_exp)
#efnb3_18 -> Nhs -> efnb3_exp
mediation_analysis(efnb3_18, "Nhs", efnb3_exp)
#efnb3_18 -> cacna1h -> efnb3_exp 
mediation_analysis(efnb3_18, "Cacna1h", efnb3_exp)
#efnb3_18 -> ghsr -> efnb3_exp
mediation_analysis(efnb3_18, "Ghsr", efnb3_exp)
#efnb3_18 -> Crhr2 -> efnb3_exp
mediation_analysis(efnb3_18, "Crhr2", efnb3_exp)
#ghsr_18 -> Nhs -> ghsr_exp
mediation_analysis(ghsr_18, "Nhs", ghsr_exp)
#ghsr_18 -> efnb3 -> ghsr_exp
mediation_analysis(ghsr_18, "Efnb3", ghsr_exp)
#ghsr_18 -> Cacna1h -> ghsr_exp
mediation_analysis(ghsr_18, "Cacna1h", ghsr_exp)
#ghsr_18 -> Kcnd2-> ghsr_exp
mediation_analysis(ghsr_18, "Kcnd2", ghsr_exp)
#ghsr_4 -> arg1 -> ghsr_exp
mediation_analysis(ghsr_18, "Arg1", ghsr_exp)

# None of these genes are satisfying the conditions for mediation analysis >:(
# With a slightly more lenient confidence level for step 3, Arg1 was a mediator for ghsr4,
# Cacna1h was a mediator for ghsr18, ghsr was a mediator for efnb318, cacna1h

anova(lm(ghsr_exp~sex + efnb3_exp + ghsr_18))


cacna1h_6 <- get_genoprob(6, 4.323241)
cacna1h_18 <- get_genoprob(18, 4.602510)
slc16a7_6 <- get_genoprob(6, 4.776124)

arg1_exp <- get_exp_dat("Arg1")
ptprz1_exp <- get_exp_dat("Ptprz1")
kcnd2_exp <- get_exp_dat("Kcnd2")
cacna1h_exp <- get_exp_dat("Cacna1h")
armc4_exp <- get_exp_dat("Armc4")
dscam_exp <- get_exp_dat("Dscam")
spock3_exp <- get_exp_dat("Spock3")
nhs_exp <- get_exp_dat("Nhs")
crhr2_exp <- get_exp_dat("Crhr2")
hhex_exp <- get_exp_dat("Hhex")
sst_exp <- get_exp_dat("Sst")
slc16a7_exp <- get_exp_dat("Slc16a7")
efnb3_exp <- get_exp_dat("Efnb3")
ly6h_exp <- get_exp_dat("Ly6h")
wnt10a_exp <- get_exp_dat("Wnt10a")
f5_exp <- get_exp_dat("F5")
ncam2_exp <- get_exp_dat("Ncam2")

# At: steps 1,2,4 = 0.01, step 3 = 0.05
mediation_analysis(cacna1h_6, ptprz1_exp, cacna1h_exp)
# Passes mediation analysis! :))
mediation_analysis(cacna1h_6, spock3_exp, cacna1h_exp)
# Does not pass
mediation_analysis(cacna1h_6, dscam_exp, cacna1h_exp)
# Does not pass
mediation_analysis(cacna1h_6, slc16a7_exp, cacna1h_exp)
# Passes
mediation_analysis(cacna1h_6, hhex_exp, cacna1h_exp)
# Does not pass
mediation_analysis(cacna1h_6, sst_exp, cacna1h_exp)
# DOES NOT PASS

mediation_analysis(cacna1h_18, armc4_exp, cacna1h_exp)
# Passes!
mediation_analysis(cacna1h_18, nhs_exp, cacna1h_exp)
# does not pass
mediation_analysis(cacna1h_18, efnb3_exp, cacna1h_exp)
# does not pass
mediation_analysis(cacna1h_18, ghrelin_list$ghsr_exp, cacna1h_exp)
# does not pass
mediation_analysis(cacna1h_18, crhr2_exp, cacna1h_exp)
# does not pass

mediation_analysis(slc16a7_6, ptprz1_exp, slc16a7_exp)
# Passed!!! :)


slc16a7_exp <- get_exp_dat("Slc16a7")
ptprz1_exp
slc16a7_6

sst_6 <- get_genoprob(6, 4.323241)
sst_18 <- get_genoprob(18, 8.263431)

mediation_analysis(sst_6, ptprz1_exp, sst_exp)
# Passes
mediation_analysis(sst_6, spock3_exp, sst_exp)
# does not pass
mediation_analysis(sst_6, dscam_exp, sst_exp)
# does not pass
mediation_analysis(sst_6, cacna1h_exp, sst_exp)
# does not pass
mediation_analysis(sst_6, hhex_exp, sst_exp)
# Passes


mediation_analysis(sst_18, cacna1h_exp, sst_exp)
# does not pass
mediation_analysis(sst_18, armc4_exp, sst_exp)
# does not pass
mediation_analysis(sst_18, ly6h_exp, sst_exp)
# does not pass
mediation_analysis(sst_18, crhr2_exp, sst_exp)
# does not pass
mediation_analysis(sst_18, wnt10a_exp, sst_exp)
# does not pass

ptprz1_6 <- get_genoprob(9, 65.580696)

mediation_analysis(ptprz1_6, f5_exp, ptprz1_exp)
# does not pass

spock3_13 <- get_genoprob(13, 98.487897)
# This wont work :((...not sure why
spock3_6 <- get_genoprob(6, 4.323241)

mediation_analysis(spock3_13, ptprz1_exp, spock3_exp)


mediation_analysis(spock3_6, ptprz1_exp, spock3_exp)
# Passes



num_expressions("Armc4")

nhs_18_1 <- get_genoprob(18, 4.644762)
nhs_18_2 <- get_genoprob(18, 43.103868)
fem1c_exp <- get_exp_dat("Fem1c")

mediation_analysis(nhs_18_1, ghrelin_list$ghsr_exp, nhs_exp)
# does not pass
mediation_analysis(nhs_18_1, cacna1h_exp, nhs_exp)
# does not pass
mediation_analysis(nhs_18_1, ghrelin_list$efnb3_exp, nhs_exp)
# does not pass

mediation_analysis(nhs_18_2, fem1c_exp, nhs_exp)
# does not pass



hhex_6 <- get_genoprob(6, 4.927084)
hhex_18 <- get_genoprob(18, 46.223048)
hhex_19 <- get_genoprob(19, 37.387306)
ghsr_exp <- ghrelin_list$ghsr_exp
rbp4_exp <- get_exp_dat("Rbp4")

mediation_analysis(hhex_6, ptprz1_exp, hhex_exp)
# passes
mediation_analysis(hhex_6, ghsr_exp, hhex_exp)
# does not pass
mediation_analysis(hhex_6, spock3_exp, hhex_exp)
# does not pass
mediation_analysis(hhex_6, dscam_exp, hhex_exp)
# does not pass
mediation_analysis(hhex_6, rbp4_exp, hhex_exp)
# PASSES

mediation_analysis(hhex_18, ptprz1_exp, hhex_exp)
# does not pass
mediation_analysis(hhex_18, nhs_exp, hhex_exp)
# does not pass
mediation_analysis(hhex_18, ncam2_exp, hhex_exp)
# does not pass

mediation_analysis(hhex_19, rbp4_exp, hhex_exp)
# Does not pass

dscam_6 <- get_genoprob(6, 4.323241)
dscam_18 <- get_genoprob(18, 7.761496)

mediation_analysis(dscam_6, ptprz1_exp, dscam_exp)
# passes
mediation_analysis(dscam_6, slc16a7_exp, dscam_exp)
# passes (at lower significance)
mediation_analysis(dscam_6, spock3_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_6, sst_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_6, cacna1h_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_6, hhex_exp, dscam_exp)
# does not pass


mediation_analysis(dscam_18, nhs_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_18, ghsr_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_18, armc4_exp, dscam_exp)
# does not pass
mediation_analysis(dscam_18, cacna1h_exp, dscam_exp)
# does not pass


cor(rbp4_exp, ghrelin_list$food_ave)

rbp4_6 <- get_genoprob(6, 13.081484)
mediation_analysis(rbp4_6, ptprz1_exp, rbp4_exp)

num_expressions("Ptprz1")

rbm12b2_exp <- get_exp_dat("Rbm12b2")
ghsr_4 <- get_genoprob(4, 11.984686)
cacna1h_4 <- get_genoprob(4, 4.570163)
efnb3_4 <- get_genoprob(4, 7.127355)


mediation_analysis(ghsr_4, rbm12b2_exp, ghsr_exp)
mediation_analysis(ghsr_4, arg1_exp, ghsr_exp)
mediation_analysis(cacna1h_4, rbm12b2_exp, cacna1h_exp)
mediation_analysis(cacna1h_4, arg1_exp, efnb3_exp)
mediation_analysis(efnb3_4, rbm12b2_exp, efnb3_exp)
mediation_analysis(efnb3_4, arg1_exp, efnb3_exp)
# NONE OF THESE pass mediation analysis




anova(lm(ghsr_exp ~ sex + ghsr_4))[2,5]
# middle is linked to begin
anova(lm(arg1_exp ~ sex + ghsr_4))[2,5]
# end NOT linked to begin after accounting for middle
anova(lm(ghsr_exp ~ sex + arg1_exp + ghsr_4))[3,5]
# middle and begin linked after accounting for end
anova(lm(rbm12b2_exp ~ sex + arg1_exp + ghsr_4))[3,5]


