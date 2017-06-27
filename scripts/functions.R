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

# function to calculate p-values for one phenotype and many variates (gene expressions
# in this case); 
# In my usages, express_dat refers to the rankz.mrna object and covar refers to sex
pvals_changing_covar <- function(pheno, express_dat, covar) {
  pvals <- numeric(0)
  for(i in 1:ncol(express_dat)) {
    pvals <- as.numeric(as.character(c(pvals, anova(lm(pheno ~ covar + express_dat[,i]))[2,5])))
    print(i)
  }
  return(p.adjust(pvals, method = "BH")) #return list of adjusted pvals 
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

# Takes a list of (named) p-values and returns a list of all under the significance level
sig_list <- function(pval_vec, siglev) {
  empty_vec <- numeric(0)
  for(i in 1:length(pval_vec)) {
    if(pval_vec[i] < siglev ) {
      empty_vec <- c(empty_vec, pval_vec[i])
    }
  }
  sorted_sig <- sort(empty_vec, decreasing = FALSE)
  return(sorted_sig)
}

# Function that takes a chromosome and position and returns
# the genoprobs at that position
get_genoprob <- function(chr, pos){ 
  chr_data <- snps[snps$chr == chr,]
  chr_data$cM <- chr_data$cM - pos
  chr_data <- chr_data[chr_data$cM >= 0,]
  location <- which(snps$marker == chr_data$marker[1])
  return(genoprobs[,,location])
}

get_gene_names <- function(){
  gene_names <- character(0)
  for(i in 1:ncol(rankz.mrna)) { 
    gene_name <- annot.mrna$symbol[i]
    gene_names <- c(gene_names, gene_name)
  }
  return(gene_names)
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

# BIC Score model analysis
# X: gene expression data for a given gene
# Y: quantitative measurement of clinical phenotype
# Q: genotype at a given marker
triple.fit <- function(X, Y, Q) {
  
  # Remove any NA values from the data
  #indx <- sort(unique(c(which(is.na(X)), which(is.na(Y)), which(is.na(Q)))))
  #X <- X[-indx]
  #Y <- Y[-indx]
  #Q <- Q[-indx]
  #print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  # Calculate BIC scores for models
  bic.independent <- BIC(lm(X~Q)) + BIC(lm(Y~Q)) # X<-Q->Y
  bic.reactive <- BIC(lm(X~Y)) + BIC(lm(Y~Q)) # Q->Y->X
  bic.causal <- BIC(lm(X~Q)) + BIC(lm(Y~X)) # Q->X->Y
  bic.complex <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))
  
  # Print out the scores from each model
  print("BIC Scores of each model")
  scores <- c(bic.independent, bic.reactive, bic.causal, bic.complex)
  names(scores) <- c("independent", "reactive", "causal", "complex")
  print(scores)
  
  # Make lowest BIC score 0 and linearize all other scores accordingly to calculate Delta values
  deltas <- scores - min(scores)
  
  # Estimate the strength of evidence for each model
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  # Print out the probabilities of each model being the likely explanation for the data
  print("Probability of each model explaining the data")
  print(strengths * 100)
  
  # Print out how many more times likely the best model is
  print("The factor by which the best model is better than the rest")
  print(max(strengths) / strengths)
  
}





























