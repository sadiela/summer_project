# Running through Gary's Script
# Sadie Allen
# July 31, 2017

# libraries
library(tidyverse)
library(corrplot)
library(ggtern)
library(cowplot)
library(GGally)

# cd and load data
setwd("/Users/s-allens/Documents/ssp/summer_project/")
load("/Users/s-allens/Documents/ssp/summer_project/cell_type_markers.RData")

###
# create a tidy dataframe the for data
expr.data <- cbind(phenotypes, alpha_markers, delta_markers, beta_markers) %>%
  as_tibble() %>%
  rename(Sex=sex, Gen=generation, BW16=weight_16wk, Days=diet_days, FoodAvg=food_ave)
# columns 6-13 alpha
# columns 14-21 delta
# columns 22-29 beta  #24-25 are Ins

###
# clean up
rm(phenotypes, alpha_markers, delta_markers, beta_markers)

###
# function to compute residuals to adjust x wrt one covariate
my.adjust1 <- function(x, cov){
  res <- residuals(lm(x ~ cov, na.action=na.exclude))
  mu <- mean(x, na.rm=TRUE)
  new_x <- mu + res
  new_x
}

###
# function to compute residuals to adjust x wrt two covariates
my.adjust2 <- function(x, cov1, cov2){
  res <- residuals(lm(x ~ cov1 + cov2, na.action=na.exclude))
  mu <- mean(x, na.rm=TRUE)
  new_x <- mu + res
  new_x
}

# log transform and adjust GEX data for Gen
expr.adj <- expr.data
for(i in 6:29){
  expr.adj[i] <- my.adjust1(log(expr.data[[i]]), expr.data$Gen)
}

###
# look at some scatter plots of GEX data
quartz()
ggplot(expr.adj, aes(x=Gcg, y=Ins1, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)

###
# plot beta cell markers against BW in one panel
quartz()
expr.adj %>% select(Sex, BW16,Pdx1:Slc2a2) %>%
  gather(key="Gene", value="Expr", Pdx1:Slc2a2) %>%
  ggplot(aes(x=Expr, y=BW16, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~Gene, nrow=2, scales="free_x")
#expression pattern is consistent across all genes except Ins1, Ins2

###
# plot alpha cell markers against BW in one panel
quartz()
expr.adj %>% select(Sex, BW16,Mafb:Sgce) %>%
  gather(key="Gene", value="Expr", Mafb:Sgce) %>%
  ggplot(aes(x=Expr, y=BW16, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~Gene, nrow=2, scales="free_x")
#expression pattern is consistent across all genes

###
# plot delta cell markers against BW in one panel
quartz()
expr.adj %>% select(Sex, BW16,Rbp4:F5) %>%
  gather(key="Gene", value="Expr", Rbp4:F5) %>%
  ggplot(aes(x=Expr, y=BW16, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~Gene, nrow=2, scales="free_x")
#expression pattern is consistent across all genes

###
# look at correlation by sex
quartz()
expr.adj %>%
  filter(Sex=="F") %>% select(4:29) %>%
  cor() %>% corrplot(title="Female")
#
quartz()
expr.adj %>%
  filter(Sex=="M") %>% select(4:29) %>%
  cor() %>% corrplot(title="Male")
#

###
# remove Sex effects before computing the PC
# note overwrite of old expr.adj
for(i in 6:29){
  expr.adj[i] <- my.adjust2(log(expr.data[[i]]), expr.data$Gen, expr.data$Sex)
}

###
# compute principle components for alpha cells
pc_alpha <- princomp(expr.adj[,6:13], cor=TRUE)
# confirm that pc1 is dominant
quartz()
plot(pc_alpha)
# check loadings
pc_alpha$loadings[,1:4]
# save pc1
expr.adj <- mutate(expr.adj, Alpha = -pc_alpha$scores[,1])

###
# compute principle components for delta cells
pc_delta <- princomp(expr.adj[,14:21], cor=TRUE)
# confirm that pc1 is dominant
quartz()
plot(pc_delta)
# check loadings
pc_delta$loadings[,1:4]
# save pc1
expr.adj <- mutate(expr.adj, Delta = -pc_delta$scores[,1])

###
# compute principle components for beta cells
pc_beta <- princomp(expr.adj[,c(22:23,26:29)], cor=TRUE)
# confirm that pc1 is dominant
quartz()
plot(pc_beta)
# not as tight
# check the loadings
pc_beta$loadings[,1:4]
# Ucn3 and Slc2a2 are mild outliers
# save pc1
expr.adj <- mutate(expr.adj, Beta = -pc_beta$scores[,1])


###
# compute principle components for insulin
pc_ins <- princomp(expr.adj[,24:25], cor=TRUE)
# confirm that pc1 is dominant
quartz()
plot(pc_ins)
# check the loadings
pc_ins$loadings[,1:2]
# save pc1
expr.adj <- mutate(expr.adj, Ins = -pc_ins$scores[,1])

###
# plots to check the PCs
quartz()
ggplot(expr.adj, aes(x=Irx2, y=Alpha, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)
#
quartz()
ggplot(expr.adj, aes(x=Ptprz1, y=Delta, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)
#
quartz()
ggplot(expr.adj, aes(x=Pdx1, y=Beta, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)
#
quartz()
ggplot(expr.adj, aes(x=Ins1, y=Ins, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)

###
# look at correlation of PCs by sex
quartz()
expr.adj %>%
  filter(Sex=="F") %>% select(BW16, FoodAvg, Alpha:Ins) %>%
  cor() %>% corrplot(title="Female")
#
quartz()
expr.adj %>%
  filter(Sex=="M") %>% select(BW16, FoodAvg, Alpha:Ins) %>%
  cor() %>% corrplot(title="Male")
#

###
# plot PC traits against BW16 in one panel
quartz()
expr.adj %>% select(Sex, BW16, Alpha:Ins) %>%
  gather(key="PC", value="Expr", Alpha:Ins) %>%
  ggplot(aes(x=Expr, y=BW16, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~PC, nrow=2, scales="free_x")

###
# plot PC traits against FoodAvg in one panel
quartz()
expr.adj %>% select(Sex, FoodAvg, Alpha:Ins) %>%
  gather(key="PC", value="Expr", Alpha:Ins) %>%
  ggplot(aes(x=Expr, y=FoodAvg, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~PC, nrow=2, scales="free_x")

###
# regression on BW16
# alpha cells trump delta cells
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + Delta + Alpha, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + Alpha + Delta, data=expr.adj))

# insulin always matters
anova(lm(BW16 ~ Sex + Gen + Beta + Alpha + Delta + Ins, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + Alpha + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Ins + Beta + Alpha + Delta, data=expr.adj))

# beta always matter
anova(lm(BW16 ~ Sex + Gen + Ins + Alpha + Delta + Beta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Ins + Alpha + Beta + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Ins + Beta + Alpha + Delta, data=expr.adj))

###
# regression on FoodAvg
# alpha trumps delta
anova(lm(FoodAvg ~ Sex + Gen + Beta + Ins + Delta + Alpha, data=expr.adj))
anova(lm(FoodAvg ~ Sex + Gen + Beta + Ins + Alpha + Delta, data=expr.adj))

# insulin always matters
anova(lm(FoodAvg ~ Sex + Gen + Beta + Alpha + Delta + Ins, data=expr.adj))
anova(lm(FoodAvg ~ Sex + Gen + Beta + Ins + Alpha + Delta, data=expr.adj))
anova(lm(FoodAvg ~ Sex + Gen + Ins + Beta + Alpha + Delta, data=expr.adj))

# alpha cells trump beta
anova(lm(FoodAvg ~ Sex + Gen + Ins + Alpha + Delta + Beta, data=expr.adj))
anova(lm(FoodAvg ~ Sex + Gen + Ins + Alpha + Beta + Delta, data=expr.adj))
anova(lm(FoodAvg ~ Sex + Gen + Ins + Beta + Alpha + Delta, data=expr.adj))

###
# regression on BW16 with FoodAvg as predictor

# FoodAvg always matters
anova(lm(BW16 ~ Sex + Gen + FoodAvg + Beta + Ins + Alpha + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + FoodAvg + Alpha + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + Alpha + Delta + FoodAvg, data=expr.adj))

# beta always matters
anova(lm(BW16 ~ Sex + Gen + Beta + Ins + FoodAvg + Alpha + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Ins + FoodAvg + Beta + Alpha + Delta, data=expr.adj))
anova(lm(BW16 ~ Sex + Gen + Ins + FoodAvg + Alpha + Delta + Beta, data=expr.adj))

#
quartz()
ggplot(expr.adj, aes(x=FoodAvg, y=BW16, color=Sex)) +
  geom_point() + geom_smooth(method="lm", se=FALSE)

