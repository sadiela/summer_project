# Principle Component Analysis
# Sadie Allen
# July 14, 2017
# Using principle component analysis to determine the directions of maximum
# variation in phenotype data; group by different characteristics

rm(list = ls())

#Load Libraries
library(tidyverse)
library(plotly)
library(pcaMethods)
library(RColorBrewer)

# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

phenos <- read.csv("data/matched_phenos.csv")
full_phenos <- phenos

sex <- as.numeric(as.factor((annot.samples$Sex))) -1 
generation <- as.numeric(as.factor(annot.samples$Generation))


x <- vector("numeric")
for(i in 1:length(full_phenos)){
  if(is.numeric(full_phenos[,i])== TRUE && all(is.na(phenos[,i])) == FALSE){
    x <- c(x,i)
  }
}
full_phenos <- full_phenos[,x]

full_phenos <- full_phenos[, -grep("DOwave", colnames(full_phenos))]
full_phenos <- full_phenos[,-1]

full_phenos <- apply(full_phenos, 2, scale)


full_pca <- pca(full_phenos, nPcs = 141)

hist(phenos$num_islets)
#pch -> plotting character
#col -> colors
col <- sex + 1
col_g <- generation
quartz()
plotPcs(full_pca, pcs = 1:4, type = c("scores", "loadings"), sl = NULL, hotelling = NULL, col = col, pch = 16)

plotPcs(full_pca, pcs = 1:4, type = c("scores", "loadings"), sl = NULL, hotelling = NULL, col = col_g, pch = 16)


quartz()

breaks = seq(min(phenos$num_islets), max(phenos$num_islets), 20)
col = colorRampPalette(brewer.pal(11, "Spectral"))(length(breaks) - 1)
plotPcs(full_pca, pcs = 1:5, type = c("scores", "loadings"), sl = NULL, hotelling = NULL, col = col, pch = 16)


biplot(full_pca, choices = 1:2, scale=1, pc.biplot = FALSE)

dim(full_pca)


plot.pcaRes(all_pca)

# Okayyyy that was funn... but I didn't reall learn anything from it because I couldn't color the
# data points


















