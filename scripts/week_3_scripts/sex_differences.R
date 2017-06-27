# Sex Differences
# Sadie Allen
# June 23, 2017
# In this script I will analyze how sex affects different phenotypes
# and gene expressions. 

## Load Libraries ##
library(ggplot2)

## Load in data ##
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)
rownames(matched_phenotypes) <- ghrelin_shortlist[,1]
ghrelin_shortlist$X <- NULL

fat_wt <- matched_phenotypes$fat_pad_weight
gcg_content <- matched_phenotypes$Gcg_content
ghrelin_shortlist <- cbind(ghrelin_shortlist, leptin, fat_wt, gcg_content)

## Subsetting data by gender ##
gender_sep <- split(ghrelin_shortlist, ghrelin_shortlist$sex)
females <- gender_sep$F
males <- gender_sep$M

## Create histograms split by gender and performing t-tests ##
# For t-tests: when specifying alternative hypothesis, "greater" means 
# first input is greater than the second

# Differences in Glu_0min 
ggplot(ghrelin_shortlist, aes(x = Glu_0min, group = sex, fill = sex)) + 
      geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$Glu_0min, males$Glu_0min, alternative = "l")
# p ≈ 0, female glu_0min < male glu_0min

ggplot(ghrelin_shortlist, aes(x = Glu_sac, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$Glu_sac, males$Glu_sac, alternative = "l")
# p ≈ 0, female glu_sac < male glu_sac

ggplot(ghrelin_shortlist, aes(x = Ins_0min, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$Ins_0min, males$Ins_0min, alternative = "l")
# NOT NORMAL.... but p ≈ 0, female ins_0min < male ins_0min

ggplot(ghrelin_shortlist, aes(x = Ins_sac, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$Ins_sac, males$Ins_sac, alternative = "l")
# p ≈ 0, female ins_sac < males ins_sac

ggplot(ghrelin_shortlist, aes(x = food_ave, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$food_ave, males$food_ave, alternative = "l")
# p ≈ 0, female food_ave < male food_ave

ggplot(ghrelin_shortlist, aes(x = weight_sac, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$weight_sac, males$weight_sac, alternative = "l")
# p ≈ 0, female weight_sac < male weight_sac

ggplot(ghrelin_shortlist, aes(x = G33_ins_secrete, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$G33_ins_secrete, males$G33_ins_secrete, alternative = "l")
# p ≈ 0, female g33 < male g33

ggplot(ghrelin_shortlist, aes(x = G83_ins_secrete, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$G83_ins_secrete, males$G83_ins_secrete, alternative = "l")
# p ≈ 0, female g83 < male g83

ggplot(ghrelin_shortlist, aes(x = G167_ins_secrete, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$G167_ins_secrete, males$G167_ins_secrete, alternative = "l")
# p ≈ 0, female g167 < male g167

ggplot(ghrelin_shortlist, aes(x = ghrl, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$ghrl, males$ghrl, alternative = "g")
# p ≈ 0, female ghrl > male ghrl

ggplot(ghrelin_shortlist, aes(x = ghsr, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$ghsr, males$ghsr, alternative = "g")
# p ≈ 0, female ghsr > male ghsr

ggplot(ghrelin_shortlist, aes(x = sst, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$sst, males$sst, alternative = "g")
# p ≈ 0, female sst > male sst

ggplot(ghrelin_shortlist, aes(x = ins1, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$ins1, males$ins1, alternative = "l")
# p ≈ 0, female ins1 < male ins1

ggplot(ghrelin_shortlist, aes(x = ins2, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$ins2, males$ins2, alternative = "l")
# p ≈ 0, female ins2 < male ins2

ggplot(ghrelin_shortlist, aes(x = sstr, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.1) + theme_bw()
t.test(females$sstr, males$sstr)
# p = 0.3425, sstr not significantly different between sexes

ggplot(ghrelin_shortlist, aes(x = num_islets, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$num_islets, males$num_islets)
# p = 0.06533, number of islets not significantly different between sexes

ggplot(ghrelin_shortlist, aes(x = Ins_per_islet, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$Ins_per_islet, males$Ins_per_islet)
# p = 0.2226, Insulin per islet not significantly different between sexes

ggplot(ghrelin_shortlist, aes(x = WPIC, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$WPIC, males$WPIC) #default is two-sided
# p = 0.07009, WPIC's not significantly different between sexes

ggplot(ghrelin_shortlist, aes(x = weight_change_ave, group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(females$weight_change_ave, males$weight_change_ave, alternative = "l")
# p ≈ 0, female wca < male wca

# Differences in fat_wt 
ggplot(ghrelin_shortlist, aes(x = log10(fat_wt), group = sex, fill = sex)) + 
  geom_histogram(position = "dodge", binwidth = 0.05) + theme_bw()
t.test(log10(females$fat_wt), log10(males$fat_wt), alternative = "l")
# p ≈ 0, female fat weight < male fat weight (doesnt look terribly normal even when transformed)

# Gender-affected phenotypes: weight_change_ave, glu_0min, glu_sac, ins_sac, 
# food_ave, weight_sac, g33_ins_secrete, g83_ins_secrete, g167_ins_secrete, 
# ghrl, ghsr, sst, ins1, ins2

## Scatterplots separated by gender ##

#ghrelin receptor vs food consumption vs bodyweight
ggplot(ghrelin_shortlist, aes(x = ghsr, y = food_ave, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$food_ave)
# -0.4745 correlation

ggplot(ghrelin_shortlist, aes(x = food_ave, y = weight_sac, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$food_ave, ghrelin_shortlist$weight_sac)
# 0.726 correlation

ggplot(ghrelin_shortlist, aes(x = ghsr, y = weight_sac, group = sex, fill = sex, col = sex)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
cor(ghrelin_shortlist$ghsr, ghrelin_shortlist$weight_sac)
# -0.509 correlation



