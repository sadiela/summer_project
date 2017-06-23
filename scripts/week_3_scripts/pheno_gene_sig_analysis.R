# Analyzing Phenotypical Gene Significance
# Sadie Allen
# June 22, 2017

# In this script, I explore the contents of the pheno_gene_exp_signif.csv
# file to see which genes seem to have significant effects on my phenotypes
# of interest.

# Load data table
pheno_gene_exp_signif <- read.csv("/Users/s-allens/Documents/ssp/summer_project/results/pheno_gene_exp_signif.csv")
pheno_gene_exp_signif$genes <- as.character(pheno_gene_exp_signif$genes)
pheno_gene_exp_signif$X <- NULL
gene_names <- pheno_gene_exp_signif[,1]

# Function:
sig_list <- function(pval_vec) {
  empty_vec <- numeric(0)
  for(i in 1:length(pval_vec)) {
    if(pval_vec[i] < 0.001 ) {
      empty_vec <- c(empty_vec, pval_vec[i])
    }
  }
  return(empty_vec)
}


# First Phenotype: Glu_0min (plasma glucose at time 0 for the oGTT collected after a 4 hr fast)
Glu_0min_pvals <- pheno_gene_exp_signif[,2]
names(Glu_0min_pvals) <- gene_names
Glu_0min_pvals <- sort(Glu_0min_pvals, decreasing = FALSE)
Glu_0min_sig_genes <- sig_list(Glu_0min_pvals)
# This gave a list of 640 significant genes! Yikes!

# Glu_sac
glu_sac_pvals <- pheno_gene_exp_signif[,3]
names(glu_sac_pvals) <- gene_names
glu_sac_pvals <- sort(glu_sac_pvals, decreasing = FALSE)
glu_sac_sig_genes <- sig_list(glu_sac_pvals)

# Ins_0min
ins_0min_pvals <- pheno_gene_exp_signif[,4]
names(ins_0min_pvals) <- gene_names
ins_0min_pvals <- sort(ins_0min_pvals, decreasing = FALSE)
ins_0min_sig_genes <- sig_list(ins_0min_pvals)

# Ins_sac
ins_sac_pvals <- pheno_gene_exp_signif[,5]
names(ins_sac_pvals) <- gene_names
ins_sac_pvals <- sort(ins_sac_pvals, decreasing = FALSE)
ins_sac_sig_genes <- sig_list(ins_sac_pvals)

# Food_ave
food_ave_pvals <- pheno_gene_exp_signif[,6]
names(food_ave_pvals) <- gene_names
food_ave_pvals <- sort(food_ave_pvals, decreasing = FALSE)
food_ave_sig_genes <- sig_list(food_ave_pvals)

# Weight_sac
weight_sac_pvals <- pheno_gene_exp_signif[,7]
names(weight_sac_pvals) <- gene_names
weight_sac_pvals <- sort(weight_sac_pvals, decreasing = FALSE)
weight_sac_sig_genes <- sig_list(weight_sac_pvals)

# G33_ins_secrete
g33_ins_secrete_pvals <- pheno_gene_exp_signif[,8]
names(g33_ins_secrete_pvals) <- gene_names
g33_ins_secrete_pvals <- sort(g33_ins_secrete_pvals, decreasing = FALSE)
g33_ins_secrete_sig_genes <- sig_list(g33_ins_secrete_pvals)

# Num_islets
num_islet_pvals <- pheno_gene_exp_signif[,9]
names(num_islet_pvals) <- gene_names
num_islet_pvals <- sort(num_islet_pvals, decreasing = FALSE)
num_islet_sig_genes <- sig_list(num_islet_pvals)

# Ins_per_islet
ins_per_islet_pvals <- pheno_gene_exp_signif[,10]
names(ins_per_islet_pvals) <- gene_names
ins_per_islet_pvals <- sort(ins_per_islet_pvals, decreasing = FALSE)
ins_per_islet_sig_genes <- sig_list(ins_per_islet_pvals)

# WPIC
wpic_pvals <- pheno_gene_exp_signif[,11]
names(wpic_pvals) <- gene_names
wpic_pvals <- sort(wpic_pvals, decreasing = FALSE)
wpic_sig_genes <- sig_list(wpic_pvals)

# Weight_change_ave
weight_change_ave_pvals <- pheno_gene_exp_signif[,12]
names(weight_change_ave_pvals) <- gene_names
weight_change_ave_pvals <- sort(weight_change_ave_pvals, decreasing = FALSE)
weight_change_ave_sig_genes <- sig_list(weight_change_ave_pvals)

