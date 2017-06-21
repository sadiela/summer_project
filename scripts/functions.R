#FUNCTION LIBRARY
#load qtl library
library(qtl)
library(ggplot2)

#functions:
#panel.cor and panel.hist: use with pair() to create gridded scatterplots
#gene_exp_data: get mrna expression data for different genes
#pval_vec: returns a vector of pvalues testing the significance of a covariate on a phenotype
#significance: sorts results of pval_vec
#fit:
#fit_c1:
#fit_c2:

# some useful plotting functions
# see documentation for "pairs" function 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.35)+1])
}
#
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

#function to get mrna expression data for different genes
#takes: name of gene
#returns: column containing normalized expression data from the 378 mice
get_exp_dat <- function(gene_name) {
  #use the annot.mrna dataframe to get the ID of the gene requested
  gene_id <- annot.mrna$id[annot.mrna$symbol == gene_name]
  #get list of data
  expr_data <- rankz.mrna[, colnames(rankz.mrna) == gene_id]
  return(expr_data)
}

#making a function that returns a list of pvalues after
#testing the importance of a covariate to a phenotype
#USAGE: dataset must contain the covariate in question 
pval_vec <- function(dataset, covar) {
  pvals <- numeric(0)
  name_vec <- character(0)
  for(i in names(dataset)) {
    #print(i)
    if(class(dataset[,colnames(dataset) == i]) == "numeric") {
      pvals <- c(pvals, fit(i, covar, dataset))
      name_vec <- c(name_vec, i)
      print(i)
    } else {
      print(paste(i,"is not numeric", sep = " "))
    }
  }
  names(pvals) <- name_vec
  return(pvals)
}

#list of pvals testing a variable when one covariate has already been identified
#USAGE: dataset must contain the covariate in question 
pval_vec_1co <- function(dataset, covar, var) {
  pvals <- numeric(0)
  name_vec <- character(0)
  for(i in names(dataset)) {
    #print(i)
    if(class(dataset[,colnames(dataset) == i]) == "numeric") {
      pvals <- c(pvals, fit_c1(i, var, covar, dataset))
      name_vec <- c(name_vec, i)
      print(i)
    } else {
      print(paste(i,"is not numeric", sep = " "))
    }
  }
  names(pvals) <- name_vec
  return(pvals)
}

#function that checks for insignificant p values
#Wow, R doesnt allow more than one return value. Thats nice.
#to use after pval_vec, could easily be switched to return the significant pvals
significance <- function(pvals) {
  covar_insignificant <- numeric(0)
  covar_significant <- numeric(0)
  for(i in 1:length(pvals)) {
    if(pvals[i] >= 0.05) {
      covar_insignificant <- c(covar_insignificant, pvals[i])
      print(paste(names(pvals[i]), "is not affected", sep = " "))
    } else {
      covar_significant <- c(covar_significant, pvals[i])
      print(paste(names(pvals[i]), "is affected", sep = " "))
    }
  }
  info <- list(covar_insignificant, covar_significant)
  names(info) <- c("insignificant_effect", "significant_effect")
  return(info)
}

##########################################
# P-value functions (Gary)
##########################################

# function to compute pvalues unadjusted
my.aov <- function(x,f){
  fit1 <- lm(x ~ f)
  fit0 <- lm(x ~ 1)
  anova(fit0, fit1)[2,6]
}
#THIS FUNCTION IS THE SAME AS: 
#Pass pheno and covar in as strings
fit <- function(pheno, covar, dataset) {
  fit1 <- anova(lm(dataset[,colnames(dataset)==pheno] ~ dataset[,colnames(dataset)==covar], data = dataset))
  return(fit1[1, 5])
}

# then, to compute p-values adjusted for one covariate:
# function to compute pvalues adjusted for one covariate
#ORDER FOR BOTH FOLLOWING FUNCTIONS: pheno, variable of interest, confirmed covar
my.aov.c1 <- function(x,f,c){
  fit1 <- lm(x ~ f+c)
  fit0 <- lm(x ~ c)
  anova(fit0, fit1)[2,6]
}
#OR
fit_c1 <- function(pheno, var, covar, dataset) {
  fit1 <- lm(dataset[,colnames(dataset)==pheno] ~ dataset[,colnames(dataset)==covar] + dataset[,colnames(dataset) == var], data = dataset)
  fit0 <- lm(dataset[,colnames(dataset)== pheno] ~ dataset[,colnames(dataset)==covar], data = dataset)
  pval <- anova(fit0, fit1)[2,6]
  return(pval)
}

# function to compute pvalues adjusted for two covariate
#ORDER: pheno, variable of interest, covar1, covar2
my.aov.c2 <- function(x,f,c1,c2){
  fit1 <- lm(x ~ f+c1+c2)
  fit0 <- lm(x ~ c1+c2)
  anova(fit0, fit1)[2,6]
}
#OR
fit_c2 <- function(pheno, var, covar1, covar2, dataset) { 
  fit1 <- lm(dataset[,colnames(dataset)==pheno] ~ dataset[,colnames(dataset)==covar1] + dataset[,colnames(dataset) == covar2] + dataset[,colnames(dataset)== var], data = dataset)
  fit0 <- lm(dataset[,colnames(dataset)== pheno] ~ dataset[,colnames(dataset)==covar1] + dataset[,colnames(dataset) == covar2], data = dataset)
  pval <-   anova(fit0, fit1)[2,6]
  return(pval)
}

#####################
#QTL Stuff
#####################
genotype <- function(chr, pos) {
  if (chr > 0 && chr < 21 && pos >= 0 && pos <= 100) {
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr=chr, pos=pos)])
  }
}

