#Week3
#In this script, I am testing the computational capacity of my lab computer (MLG-CCMBP01)
#setwd("/Users/s-allens/Documents/ssp/summer_project")

#load libraries
library(tidyverse)
library(plotly)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)

#load important functions
source("scripts/functions.R")

#load data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)

#refine to phenotypes that I am working with 
ghrelin_shortlist <- matched_phenotypes %>% select(sex, diet_days, Glu_0min,  Glu_sac, Ins_0min, 
                         Ins_sac, food_ave, weight_sac, G33_ins_secrete, 
                         G83_ins_secrete, G167_ins_secrete)
for(i in 3:10) {
  ghrelin_shortlist[i] <- log10(ghrelin_shortlist[,i])
}

#add expression data
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

#combine with phenotypes to get final dataset
ghrelin_shortlist <- cbind(ghrelin_shortlist, gene_expressions)

#add diet_days to annot.samples
diet_days <- ghrelin_shortlist$diet_days
annot.samples <- cbind(annot.samples, diet_days)

#NOW LETS DO QTL ANALYSIS! :))
#converts DOQTL genotype probabilities to R/qtl2 format
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

#convert a marker map organized as data frame to a list
map <- map_df_to_list(map = snps, pos_column = "bp")

rownames(matched_phenotypes) <- rownames(expr.mrna)

#creating a covariate that encompasses sex, generation, and diet days
covar2 <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)

#rownames should be in order...

fit2 <- scan1(genoprobs=probs, 
              log10(ghrelin_shortlist[,3, drop = FALSE]), 
              addcovar=covar2, 
              cores=4, 
              reml=TRUE)
plot(fit2, map = map)
#############################################################################

#Now: Test Computation Power
hundred_mice <- ghrelin_shortlist[1:100,]
two_hundred_mice <- ghrelin_shortlist[1:200,]
three_hundred_mice <- ghrelin_shortlist[1:300,]

#One Phenotype:
fit_100_one <- scan1(genoprobs=probs, 
              log10(hundred_mice[,4, drop = FALSE]), 
              addcovar=covar2, 
              cores=4, 
              reml=TRUE)
#~1.86 seconds

fit_200_one <- scan1(genoprobs=probs, 
                     log10(two_hundred_mice[,4, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~4.13 seconds

fit_300_one <- scan1(genoprobs=probs, 
                     log10(three_hundred_mice[,4, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~6.36 seconds

fit_378_one <- scan1(genoprobs=probs, 
                     log10(ghrelin_shortlist[,4, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~10.16 seconds
single_pheno_times <- c(1.86, 4.13, 6.36, 10.16)

#par(mfrow=c(2,2)) (this caused a crash)
plot(fit_100_one, map = map)
plot(fit_200_one, map = map)
plot(fit_300_one, map = map)
plot(fit_378_one, map = map)

#Two Phenotypes:
fit_100_two <- scan1(genoprobs=probs, 
                     log10(hundred_mice[,4:5, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~2.94 seconds

fit_200_two <- scan1(genoprobs=probs, 
                     log10(two_hundred_mice[,4:5, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~3.53 seconds

fit_300_two <- scan1(genoprobs=probs, 
                     log10(three_hundred_mice[,4:5, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~7.14 seconds

fit_378_two <- scan1(genoprobs=probs, 
                     log10(ghrelin_shortlist[,4:5, drop = FALSE]), 
                     addcovar=covar2, 
                     cores=4, 
                     reml=TRUE)
#~14.55 seconds
two_pheno_times <- c(2.94, 3.53, 7.14, 14.55)

#par(mfrow=c(2,2)) (this caused a crash)
plot(fit_100_two, map = map) #how do I specify which phenotype I want to see?
plot(fit_200_two, map = map)
plot(fit_300_two, map = map)
plot(fit_378_two, map = map)

#Three Phenotypes:
fit_100_three <- scan1(genoprobs=probs, 
                       log10(hundred_mice[,4:6, drop = FALSE]), 
                       addcovar=covar2, 
                       cores=4, 
                       reml=TRUE)
#~2.54 seconds

fit_200_three <- scan1(genoprobs=probs, 
                       log10(two_hundred_mice[,4:6, drop = FALSE]), 
                       addcovar=covar2, 
                       cores=4, 
                       reml=TRUE)
#~3.43 seconds

fit_300_three <- scan1(genoprobs=probs, 
                       log10(three_hundred_mice[,4:6, drop = FALSE]), 
                       addcovar=covar2, 
                       cores=4, 
                       reml=TRUE)
#~7.24 seconds

fit_378_three <- scan1(genoprobs=probs, 
                       log10(ghrelin_shortlist[,4:6, drop = FALSE]), 
                       addcovar=covar2, 
                       cores=4, 
                       reml=TRUE)
#~25.01 seconds
three_pheno_times <- c(2.54, 3.43, 7.24, 25.01)

#par(mfrow=c(2,2)) (this caused a crash)
plot(fit_100_three, map = map) #how do I specify which phenotype I want to see?
plot(fit_200_three, map = map)
plot(fit_300_three, map = map)
plot(fit_378_three, map = map)

#Four Phenotypes:
fit_100_four <- scan1(genoprobs=probs, 
                      log10(hundred_mice[,4:7, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~2.50 seconds

fit_200_four <- scan1(genoprobs=probs, 
                      log10(two_hundred_mice[,4:7, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~3.36 seconds

fit_300_four <- scan1(genoprobs=probs, 
                      log10(three_hundred_mice[,4:7, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~6.43 seconds

fit_378_four <- scan1(genoprobs=probs, 
                      log10(ghrelin_shortlist[,4:7, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~14.84
four_pheno_times <- c(2.50, 3.36, 6.43, 14.84)

#par(mfrow=c(2,2)) (this caused a crash)
plot(fit_100_four, map = map) #how do I specify which phenotype I want to see?
plot(fit_200_four, map = map)
plot(fit_300_four, map = map)
plot(fit_378_four, map = map)

#Seven Phenotypes:
fit_100_seven <- scan1(genoprobs=probs, 
                      log10(hundred_mice[,4:10, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~3.40 seconds

fit_200_seven <- scan1(genoprobs=probs, 
                      log10(two_hundred_mice[,4:10, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~4.67 seconds

fit_300_seven <- scan1(genoprobs=probs, 
                      log10(three_hundred_mice[,4:10, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~16.43 seconds

fit_378_seven <- scan1(genoprobs=probs, 
                      log10(ghrelin_shortlist[,4:10, drop = FALSE]), 
                      addcovar=covar2, 
                      cores=4, 
                      reml=TRUE)
#~59.80 seconds

seven_pheno_times <- c(3.40, 4.67, 16.43, 59.80)

#par(mfrow=c(2,2)) (this caused a crash)
plot(fit_100_seven, map = map) #how do I specify which phenotype I want to see?
plot(fit_200_seven, map = map)
plot(fit_300_seven, map = map)
plot(fit_378_seven, map = map)

#No I want to graphically represent the time data I have collected. I have:
single_pheno_times
two_pheno_times
three_pheno_times
four_pheno_times
seven_pheno_times
#want a graph with time in seconds on the y axis and number of phenotypes on the x
#number of mice will be represented by different colors
times <- c(single_pheno_times, two_pheno_times, three_pheno_times, four_pheno_times, seven_pheno_times)
number_of_phenos <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), rep(7, 4))
number_of_mice <- rep(c(100,200,300,378), 5)

#combine these columns into a dataframe
scan_times <- data.frame(runtime = times, num_pheno = number_of_phenos, num_mice = number_of_mice)

#convert num_mice to a factor variable
scan_times$num_mice <- factor(scan_times$num_mice)

#scatterplot!!!
ggplot(data = scan_times, aes(x=num_pheno, y=runtime, color=num_mice)) + geom_point() + geom_line()

#yayy that was fun!


#this is where i stopped and was unable to continue or understand
#the rest of the script
# DO NOT KNOW HOW TO USE THIS
#####################################################################################
maxlod <- apply(fit2$lod, 2, max)
maxlod[maxlod>7.18]

pdf("~/phenotypes/pheno_lodplots.pdf", width=10)
for (i in setdiff(1:ncol(fit2$lod),6)) { # problems with Ins_iAUC 
  plot(fit2, lodcolumn = i, main = colnames(pheno)[i+7], ylim=c(0,11))
  abline(h=7.18, col="red", lty=2)
}  
dev.off()

png(file = "~/phenotypes/lodplot_Ins_tAUC.png", width=8.06, height = 5.71, units = "in", res=300)
plot(fit2, lodcolumn = 5, main = "Ins_tAUC")
dev.off()

coef <- scan1coef(genoprobs=probs[,11], 
                  kinship=Glist[11], 
                  pheno=log10(pheno[,12,drop=FALSE]), 
                  addcovar=covar2, 
                  reml=TRUE)
plot_coefCC(coef)

annot <- annot.mrna
annot$chr[annot$chr=="MT"] <- "M"
names(annot)[names(annot) == "middle_point"] <- "pos"

qtl.pos = which.max(fit2$lod[,5])
qtl.geno = genoprobs[,,qtl.pos]
chr = modules$chr[i]

med <- mediation.scan(modexpr[,m], 
                      mediator = expr.mrna, 
                      annot, 
                      qtl.geno = qtl.geno, 
                      covar = as.matrix(covar))