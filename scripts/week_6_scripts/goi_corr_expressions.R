# Other Genes Highly Correlated with Ghsr, Efnb3, and Svil
# Sadie Allen
# July 10, 2017
# In this script I will generate three lists, each containing genes significant
# to Ghsr, Efnb3, and Svil, respectively. The genes in each list will be ordered by
# p-value.

rm(list = ls())

# Load data
# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# Load functions
source("scripts/functions.R")

ghsr_exp <- ghrelin_list$ghsr_exp
efnb3_exp <- get_exp_dat("Efnb3")
svil_exp <- get_exp_dat("Svil")

sex <- annot.samples$Sex

gene_names <- get_gene_names()

ghsr_pvals <- pvals_changing_covar(ghsr_exp, rankz.mrna, sex)
names(ghsr_pvals) <- gene_names
ghsr_pvals <- sort(ghsr_pvals, decreasing = FALSE)
ghsr_sig_genes <- sig_list(ghsr_pvals, siglev = .0001)
head(ghsr_sig_genes)

efnb3_pvals <- pvals_changing_covar(efnb3_exp, rankz.mrna, sex)
names(efnb3_pvals) <- gene_names
efnb3_pvals <- sort(efnb3_pvals, decreasing = FALSE)
efnb3_sig_genes <- sig_list(efnb3_pvals, siglev = .0001)
head(efnb3_sig_genes)

svil_pvals <- pvals_changing_covar(svil_exp, rankz.mrna, sex)
names(svil_pvals) <- gene_names
svil_pvals <- sort(svil_pvals, decreasing = FALSE)
svil_sig_genes <- sig_list(svil_pvals, siglev = .0001)
head(svil_pvals)

save(ghsr_sig_genes, efnb3_sig_genes, svil_sig_genes, file = "data/proj_gene_exp_sig_genes.RData")
