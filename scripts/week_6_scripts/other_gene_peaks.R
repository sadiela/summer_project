# Other Gene Peaks
# Sadie Allen
# July 10, 2017
# In this script I will explore genes under some other QTL peaks: the food_ave peak on chr 7, 
# the ghsr_exp peak on chr 4, and the efnb3_exp peaks on chr 4, 11, and 18

# Utilize eQTL object

eQTLs <- read.csv(file = "data/eQTLAllAdditive.csv")

efnb3_peak_info <- eQTLs[eQTLs$symbol == "Efnb3",]
# chr 4, pos 5.566854, ci 0.0207518 - 13.32408
# chr 11, pos 19.23, ci 17.680768 - 19.61891
# chr 18, pos 3.351, ci 0.026149 - 10.25686

ghsr_peak_info <- eQTLs[eQTLs$symbol == "Ghsr",]
# Peak: chr 18, pos 4.474313, ci: 0.026149-10.21129
# Peak: chr 4, 11.984686, ci: 0.207518-13.32408

# Food_ave:
# Peak: chr 1, pos 172.2283, ci 171.3744-172.4475
# Peak: chr 7, pos 136.4545, ci = 135.5631-136.9891

# I got lists of genes for all of these peaks on ensembl.org

# Now: LOAD IN GENES
efnb3_4_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/efnb3_4_gup.csv")
efnb3_11_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/efnb3_11_gup.csv")
efnb3_18_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/efnb3_18_gup.csv")

ghsr_4_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/ghsr_4_gup.csv")

foodave_1_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/foodave_1_gup.csv")
foodave_7_genes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/all/foodave_7_gup.csv")


# chose only protein coding genes
efnb3_4_pc <- efnb3_4_genes[efnb3_4_genes$Gene.type == "protein_coding", ]
efnb3_11_pc <- efnb3_11_genes[efnb3_11_genes$Gene.type == "protein_coding", ]
efnb3_18_pc <- efnb3_18_genes[efnb3_18_genes$Gene.type == "protein_coding", ]

ghsr_4_pc <- ghsr_4_genes[ghsr_4_genes$Gene.type == "protein_coding",]

foodave_1_pc <- foodave_1_genes[foodave_1_genes$Gene.type == "protein_coding",]
foodave_7_pc <- foodave_7_genes[foodave_7_genes$Gene.type == "protein_coding",]

# Save to new data files
write.csv(efnb3_4_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/efnb3_4_pc.csv")
write.csv(efnb3_11_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/efnb3_11_pc.csv")
write.csv(efnb3_18_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/efnb3_18_pc.csv")

write.csv(ghsr_4_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/ghsr_4_pc.csv")

write.csv(foodave_1_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/foodave_1_pc.csv")
write.csv(foodave_7_pc, file = "/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/foodave_7_pc.csv")

# Efnb3
# Load list of genes w significant effects on Efnb3
load("data/proj_gene_exp_sig_genes.RData")
efnb3_sig_genes <- names(efnb3_sig_genes)

efnb3_4_pc <- efnb3_4_pc$MGI.symbol
efnb3_11_pc <- efnb3_11_pc$MGI.symbol
efnb3_18_pc <- efnb3_18_pc$MGI.symbol

intersect(efnb3_sig_genes, efnb3_4_pc) 
# Trp53inp1, Rad54b, Triqk, Plekhf2

intersect(efnb3_sig_genes, efnb3_11_pc)
# NO GENES

intersect(efnb3_sig_genes, efnb3_18_pc)
# Map3k8

# Ghsr
ghsr_sig_genes <- names(ghsr_sig_genes)

ghsr_4_genes_pc <- ghsr_4_genes_pc$MGI.symbol
ghsr_18_pc <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/genes_under_peaks/protein_coding/ghsr_18_gup_pc.csv")
ghsr_18_pc <- ghsr_18_pc$MGI.symbol

intersect(ghsr_sig_genes, ghsr_4_genes_pc)
# Rbm12b2, 1110037F02Rik, Plag1, Ints8, Rbm12b1, Ccne2, Tgs1, Clvs1, Nsmaf, Gm11808, Ubxn2b, Pdp1, Cdh17, Impad1, Asph

intersect(ghsr_sig_genes, ghsr_18_pc)
# Thoc1, Ccny, Cul2, Mtpap, Usp14, Rpl27-ps3, Fzd8, Svil

# Foodave
load("data/proj_pheno_sig_genes.RData")
food_ave_sig_genes <- names(food_ave_sig_genes)
foodave_1_genes_pc <- foodave_1_genes_pc$MGI.symbol
foodave_7_genes_pc <- foodave_7_genes_pc$MGI.symbol

intersect(food_ave_sig_genes, foodave_1_genes_pc)
# No genes

intersect(food_ave_sig_genes, foodave_7_genes_pc)
# Mki67




