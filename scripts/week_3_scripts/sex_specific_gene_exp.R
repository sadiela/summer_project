# Determining Sex-Specific Gene Expression
# Sadie Allen
# June 20, 2017

# This script is use to loop through gene expression data for over 
# 20,000 genes and determine which are significantly affected by sex


########################################################################
#                            FUNCTIONS                                 #
########################################################################

# Function that computes p-values unadjusted
# Usage: pass pheno and covar in as strings
fit <- function(pheno, covar, dataset) {
  fit1 <- anova(lm(dataset[,colnames(dataset)==pheno] ~ dataset[,colnames(dataset)==covar], data = dataset))
  return(fit1[1, 5])
}

# Function that returns a list of pvalues after
# Testing the importance of a covariate to a phenotype
# USAGE: dataset must contain the covariate in question 
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

########################################################################

# Load Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")

# Convert mrna expression data from a matrix to a dataframe
mrna_z_data <- data.frame(rankz.mrna)

# Add sex data to gene expression dataframe
sex <- annot.samples$Sex
rankz_and_gender <- cbind(sex, mrna_z_data)

# Change all variables to numerical
for(i in 2:ncol(rankz_and_gender)) {
  rankz_and_gender[,i] <- as.numeric(as.character(rankz_and_gender[,i]))
}

# Generate list of P-values
gene_exp_pvals <- pval_vec(rankz_and_gender, "sex")

# Generate list of adjusted p-values
adj_pvals <- p.adjust(gene_exp_pvals, method = "BH")

# Get the names of the genes
gene_names <- character(0)
for(i in 1:length(gene_exp_pvals)) {
  gene_name <- annot.mrna$symbol[annot.mrna$id == names(gene_exp_pvals[i])]
  gene_names <- c(gene_names, gene_name)
  names(gene_exp_pvals[i]) <- gene_name # I can't remember what this line does
}

# Create list that says if each gene is significantly affected based on its
# adjusted p-value; uneccessary, but good if you want to quickly look through
# I chose a 1% false discovery rate; 1 = sex was significant, 0 = it was not
significant <- numeric(0)
for(i in 1:length(adj_pvals)) {
  if(adj_pvals[i] < 0.01) {
    significant <- c(significant, 1)
  } else {
    significant <- c(significant, 0)
  }
}

# Build dataframe with all information
sex_dependence <- data.frame(name = gene_names, pval = gene_exp_pvals, adj_pval = adj_pvals, significance = significant)

# Sort dataframe from lowest to highest p-value
sex_dependence <- sex_dependence[order(adj_pval),]

# Save dataframe for future reference
write.csv(sex_dependence, file = "data/gene_expression_sex_dependence.csv")

# To load in the future:
sex_dependence <- read.csv("data/gene_expression_sex_dependence.csv")
