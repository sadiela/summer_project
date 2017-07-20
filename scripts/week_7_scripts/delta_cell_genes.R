# Delta Cell Genes
# Sadie Allen
# July 18, 2017
# Mark Keller provided me with a gene module that was generated based on the correlation
# structures of all the genes in the genome. This particular module includes many genes known 
# to be expressed almost exclusively in delta cells. He brought the module to my attention
# because he noticed two genes that I have been looking at, Ptprz1 and Ghsr, were included in 
# it. Upon closer inspection, almost every gene I have come across in my analysis can be found
# in the module. When Mark and Petr generated a QTL plot for the eigentrait of this module, 
# they found peaks on chromosomes 4, 6, and 18, the three peaks I have been focusing on 
# recently in my research. It is also important to note that all of the traits on the list 
# that I have analyzed thus far are negatively correlated with food intake and body weight
# phenotypes. The existence and nature of this module has several implications in my research.
# First of all, it suggests that the composition of pancreatic islets is a heritable trait
# (in the future I will run heritability estimates on the delta cell eigengene to confirm this). 
# It also broadens my hypothesis; now I am not just looking for a link between Ghsr expression 
# in the islet and food intake, I am interested in the relationship between the number of delta
# cells in the islet and food intake. 
# Here I will do more research on the genes in the module to confirm my correlation hypothesis.

#Load data
#pheno data
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# module pheno data
yellowgreen_phenos <- read.csv(file = "data/yellowgreen_phenos.csv")

# Load functions
source("scripts/functions.R")

# Top 19 genes identified in module: Ptprz1, Ghsr, Hhex, Rbp4, Gap43, F5, Nhs, Efnb3, Cacna1h
# Dscam, Cacna2d3, Fam83f, Fgf14, Slc16a7, Sst, Spock3, Mest, Arg1, Cacna1e

ptprz1_exp <- get_exp_dat("Ptprz1")
ghsr_exp <- get_exp_dat("Ghsr")
hhex_exp <- get_exp_dat("Hhex")
rbp4_exp <- get_exp_dat("Rbp4")
gap43_exp <- get_exp_dat("Gap43")
f5_exp <- get_exp_dat("F5")
nhs_exp <- get_exp_dat("Nhs")
efnb3_exp <- get_exp_dat("Efnb3")
cacna1h_exp <- get_exp_dat("Cacna1h")
dscam_exp <- get_exp_dat("Dscam")
cacna2d3_exp <- get_exp_dat("Cacna2d3")
fam83f_exp <- get_exp_dat("Fam83f")
fgf14_exp <- get_exp_dat("Fgf14")
slc16a7_exp <- get_exp_dat("Slc16a7")
sst_exp <- get_exp_dat("Sst")
spock3_exp <- get_exp_dat("Spock3")
mest_exp <- get_exp_dat("Mest")
arg1_exp <- get_exp_dat("Arg1")
cacna1e_exp <- get_exp_dat("Cacna1e")

yellowgreen_phenos <- cbind(yellowgreen_phenos, ptprz1_exp, ghsr_exp, hhex_exp, rbp4_exp, gap43_exp,
                            f5_exp, nhs_exp, efnb3_exp, cacna1h_exp, dscam_exp, cacna2d3_exp, fam83f_exp,
                            fgf14_exp, slc16a7_exp, sst_exp, spock3_exp, mest_exp, arg1_exp, cacna1e_exp)


pcor <- cor(yellowgreen_phenos, use= "complete.obs")
round(pcor, digits=2)
quartz()
corrplot(pcor, order = "hclust")































