# Significance of Gene Expressions on a List of Selected Phenotypes
# Sadie Allen
# June 21, 2017

# In this script, I will generate a dataframe containing the adjusted 
# p-values for the effects of all 20,000+ gene expressions on each phenotype
# in a list I chose based on interest and relevance to my project

########################################################################
#                            FUNCTIONS                                 #
########################################################################

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
# this function works

# In combination with pvals_changing_covar, should create a dataframe containing the 
# p-values for many phenotypes with many genes as covariates; make sure to remove identifier
# columns before using this function
# Covar = sex, express_dat = rankz.mrna, pheno_list = list of desired phenotypes, 
# gnames = list of gene names
pval_dataframe <- function(pheno_list, express_dat, covar, gnames) {
  pval_df <- data.frame(gene=gnames)
  name_vec <- c("genes")
  for(i in 1:length(pheno_list)) {
    name_vec <- c(name_vec, names(pheno_list[i]))
    pvals <- pvals_changing_covar(pheno_list[,i], express_dat, covar)
    pval_df <- cbind(pval_df, pvals)
  }
  names(pval_df) <- name_vec
  return(pval_df)
}

########################################################################

# Load Data
# mRNA expression data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")

# Phenotype data
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv", as.is=TRUE)


# I am going to have to modify my previous for loop so I can loop
# through many different covariates for one phenotype (the other script
# looped through many phenotypes with one covariate)

# Remove some phenotypes from my dataframe
mod_ghrelin_shortlist <- ghrelin_shortlist
mod_ghrelin_shortlist$sex <- NULL
mod_ghrelin_shortlist$G83_ins_secrete <- NULL
mod_ghrelin_shortlist$G167_ins_secrete <- NULL
mod_ghrelin_shortlist$X <- NULL
mod_ghrelin_shortlist$Mouse.ID <- NULL
mod_ghrelin_shortlist$ghsr <- NULL
mod_ghrelin_shortlist$sstr <- NULL
mod_ghrelin_shortlist$ins2 <- NULL
mod_ghrelin_shortlist$ghrl <- NULL
mod_ghrelin_shortlist$ins1 <- NULL
mod_ghrelin_shortlist$sst <- NULL
mod_ghrelin_shortlist$diet_days <- NULL

# adding sex to the phenotype dataframe 
sex <- annot.samples$Sex

# This will be the first column in the dataframe: a list of the gene names
gene_names <- character(0)
for(i in 1:ncol(rankz.mrna)) { 
  gene_name <- annot.mrna$symbol[i]
  gene_names <- c(gene_names, gene_name)
}

# Create the dataframe! :)
pvals_data <- pval_dataframe(mod_ghrelin_shortlist, rankz.mrna, sex, gene_names)

# Now, I will save this data frame as a csv file to be referenced in the
# future. 
write.csv(pvals_data, file = "/Users/s-allens/Documents/ssp/summer_project/results/pheno_gene_exp_signif.csv")
# To load this in the future:
pheno_gene_exp_signif <- read.csv("/Users/s-allens/Documents/ssp/summer_project/results/pheno_gene_exp_signif.csv")






