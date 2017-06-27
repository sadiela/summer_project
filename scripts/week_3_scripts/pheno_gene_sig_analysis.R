# Analyzing Phenotypical Gene Significance
# Sadie Allen
# June 22, 2017

# In this script, I explore the contents of the pheno_gene_exp_signif.csv
# file to see which genes seem to have significant effects on my phenotypes
# of interest.

load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# phenotype data
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
rownames(matched_phenotypes) <- ghrelin_shortlist[,1]


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
head(Glu_0min_sig_genes)
# Dapl1, Fkbp11, Sept6, Cngb3, Fam114a1, Nucb2
Dapl1_exp <- get_exp_dat("Dapl1")
Fkbp11_exp <- get_exp_dat("Fkbp11")
Sept6_exp <- get_exp_dat("Sept6")
Cngb3_exp <- get_exp_dat("Cngb3")
Fam114a1_exp <- get_exp_dat("Fam114a1")
Nucb2_exp <- get_exp_dat("Nucb2")

ghrelin_shortlist <- cbind(ghrelin_shortlist, Dapl1_exp, Fkbp11_exp, Sept6_exp, Cngb3_exp,
                           Fam114a1_exp, Nucb2_exp)
ggplot(ghrelin_shortlist, aes(x = Dapl1_exp, y = Glu_0min)) + geom_point()
cor(ghrelin_shortlist$Glu_0min, ghrelin_shortlist$Dapl1_exp)
#r = 0.5186

ggplot(ghrelin_shortlist, aes(x = Fkbp11_exp, y = Glu_0min)) + geom_point()
cor(ghrelin_shortlist$Glu_0min, ghrelin_shortlist$Fkbp11_exp)
#r = 0.4523

ggplot(ghrelin_shortlist, aes(x = Sept6_exp, y = Glu_0min)) + geom_point()
cor(ghrelin_shortlist$Glu_0min, ghrelin_shortlist$Sept6_exp)
#r = -0.4701

ggplot(ghrelin_shortlist, aes(x = Cngb3_exp, y = Glu_0min)) + geom_point()
cor(ghrelin_shortlist$Glu_0min, ghrelin_shortlist$Cngb3_exp)
#r = 0.4182

ggplot(ghrelin_shortlist, aes(x = Fam114a1_exp, y = Glu_0min)) + geom_point()
cor(ghrelin_shortlist$Glu_0min, ghrelin_shortlist$Fam114a1_exp)
#r = 0.4274

# Glu_sac
glu_sac_pvals <- pheno_gene_exp_signif[,3]
names(glu_sac_pvals) <- gene_names
glu_sac_pvals <- sort(glu_sac_pvals, decreasing = FALSE)
glu_sac_sig_genes <- sig_list(glu_sac_pvals)
head(glu_sac_sig_genes)
#Fkbp11, Nbr1, Osbpl2, Slc17a9, Ppapdc1b, Pex11g
Fkbp11_exp <- get_exp_dat("Fkbp11") #already there from the first pheno
Nbr1_exp <- get_exp_dat("Nbr1")
Osbpl2_exp <- get_exp_dat("Osbpl2")
Slc17a9_exp <- get_exp_dat("Slc17a9")
Ppapdc1b_exp <- get_exp_dat("Ppapdc1b")
Pex11g_exp <- get_exp_dat("Pex11g")

ghrelin_shortlist <- cbind(ghrelin_shortlist, Nbr1_exp, Osbpl2_exp, 
                           Slc17a9_exp, Ppapdc1b_exp, Pex11g_exp)

ggplot(ghrelin_shortlist, aes(x = Fkbp11_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Fkbp11_exp)
#r = 0.4577

ggplot(ghrelin_shortlist, aes(x = Nbr1_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Nbr1_exp)
#r = 0.4326

ggplot(ghrelin_shortlist, aes(x = Osbpl2_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Osbpl2_exp)
#r = 0.4530

ggplot(ghrelin_shortlist, aes(x = Slc17a9_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Slc17a9_exp)
#r = 0.4319

ggplot(ghrelin_shortlist, aes(x = Ppapdc1b_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Ppapdc1b_exp)
#r = 0.4472

ggplot(ghrelin_shortlist, aes(x = Pex11g_exp, y = Glu_sac)) + geom_point()
cor(ghrelin_shortlist$Glu_sac, ghrelin_shortlist$Pex11g_exp)
#r = 0.4124

# Ins_0min
ins_0min_pvals <- pheno_gene_exp_signif[,4]
names(ins_0min_pvals) <- gene_names
ins_0min_pvals <- sort(ins_0min_pvals, decreasing = FALSE)
ins_0min_sig_genes <- sig_list(ins_0min_pvals)
head(ins_0min_sig_genes)
# Fkbp11, Rrm2, Nucb2, Klf11, Klhl24
Fkbp11_exp <- get_exp_dat("Fkbp11") #already there from the first pheno
Rrm2_exp <- get_exp_dat("Rrm2")
Nucb2_exp <- get_exp_dat("Nucb2") # already there from a previous pheno
Klf11 <- get_exp_dat("Klf11")
Klhl24_exp <- get_exp_dat("Klhl24")
Ssr4_exp <- get_exp_dat("Ssr4")

no_0_ins2 <- ghrelin_shortlist[ghrelin_shortlist$Ins_0min > -1.3979,]

ggplot(no_0_ins2, aes(x = Fkbp11_exp, y = Ins_0min)) + geom_point()
cor(no_0_ins2$Ins_0min, no_0_ins2$Fkbp11_exp)
#r = 0.5720894

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

