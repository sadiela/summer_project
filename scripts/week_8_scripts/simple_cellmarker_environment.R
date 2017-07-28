# R Environment for Gary
# Sadie Allen
# July 27, 2017
# Gary requested a simple r environment containing:
# three data frames: one for each cell type (alpha, beta, and delta) containing expression data (RAW, NOT RANK.Z) for eight marker genes each 
# another data frame containing information about sex, generation, and diet days
# I will save this in an RData object and send it to him by the end of the day

rm(list = ls())

# Load data
# islet gene expression data, annots
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# phenotype data (need for diet days)
matched_phenos <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
rownames(matched_phenos) <- matched_phenos$Mouse.ID

#################################
# FUNCTIONS
#################################
# function to get raw mrna expression for genes: 
get_raw_exp_dat <- function(gene_name) {
  #use the annot.mrna dataframe to get the ID of the gene requested
  gene_id <- annot.mrna$id[annot.mrna$symbol == gene_name]
  #get list of data
  expr_data <- raw.mrna[, colnames(raw.mrna) == gene_id]
  return(expr_data)
}

# Function that prints the minimum and the average number of transcripts
# for the expression of a gene
num_expressions <- function(gene_name) {
  gene_id <- annot.mrna$id[annot.mrna$symbol == gene_name]
  expr_data <- raw.mrna[, colnames(raw.mrna) == gene_id]
  info <- c(min(expr_data), mean(expr_data))
  names(info) <- c("minimum expression", "average expression")
  return(info)
}

# Beta Cell Markers: Pdx-1, Glp1r, Ins1, Ins2, Mafa, Nkx6-1, Ucn3, Slc2a2 (note pdx1 is also present in delta cells)
Pdx1_exp <- get_raw_exp_dat("Pdx1")
Glp1r_exp <- get_raw_exp_dat("Glp1r")
Ins1_exp <- get_raw_exp_dat("Ins1")
Ins2_exp <- get_raw_exp_dat("Ins2")
Mafa_exp <- get_raw_exp_dat("Mafa")
Nkx61_exp <- get_raw_exp_dat("Nkx6-1")
Ucn3_exp <- get_raw_exp_dat("Ucn3")
Slc2a2_exp <- get_raw_exp_dat("Slc2a2")

beta_markers <- data.frame(Pdx1 = Pdx1_exp, Glp1r = Glp1r_exp, Ins1 = Ins1_exp, Ins2 = Ins2_exp, Mafa = Mafa_exp, Nkx61 = Nkx61_exp,
       Ucn3 = Ucn3_exp, Slc2a2 = Slc2a2_exp)
rownames(beta_markers) <- annot.samples$Mouse.ID

# Alpha Cell Markers: Mafb, Arx, Irx2, Sstr2, Mamld1, Galnt13, Gcg, Sgce
Mafb_exp <- get_raw_exp_dat("Mafb")
Arx_exp <- get_raw_exp_dat("Arx")
Irx2_exp <- get_raw_exp_dat("Irx2")
Sstr2_exp <- get_raw_exp_dat("Sstr2")
Mamld1_exp <- get_raw_exp_dat("Mamld1")
Galnt13_exp <- get_raw_exp_dat("Galnt13")
Gcg_exp <- get_raw_exp_dat("Gcg")
Sgce_exp <- get_raw_exp_dat("Sgce")

alpha_markers <- data.frame(Mafb = Mafb_exp, Arx = Arx_exp, Irx2 = Irx2_exp, Sstr2 = Sstr2_exp, Mamld1 = Mamld1_exp,
                            Galnt13 = Galnt13_exp, Gcg = Gcg_exp, Sgce = Sgce_exp)
rownames(alpha_markers) <- annot.samples$Mouse.ID


# Delta Cell Markers: Rbp4, Ptprz1, Hhex, Ghsr, Nhs, Efnb3, Gap43, F5
Rbp4_exp <- get_raw_exp_dat("Rbp4")
Ptprz1_exp <- get_raw_exp_dat("Ptprz1")
Hhex_exp <- get_raw_exp_dat("Hhex")
Ghsr_exp <- get_raw_exp_dat("Ghsr")
Nhs_exp <- get_raw_exp_dat("Nhs")
Efnb3_exp <- get_raw_exp_dat("Efnb3")
Gap43_exp <- get_raw_exp_dat("Gap43")
F5_exp <- get_raw_exp_dat("F5")

delta_markers <- data.frame(Rbp4 = Rbp4_exp, Ptprz1 = Ptprz1_exp, Hhex = Hhex_exp, Ghsr = Ghsr_exp, Nhs = Nhs_exp,
                            Efnb3 = Efnb3_exp, Gap43 = Gap43_exp, F5 = F5_exp)
rownames(delta_markers) <- annot.samples$Mouse.ID

# Sex, Generation, Diet Days object
sex_gen_diet <- data.frame(sex = annot.samples$Sex, generation = annot.samples$Generation, diet_days = matched_phenos$diet_days)
rownames(sex_gen_diet) <- annot.samples$Mouse.ID

save(alpha_markers, beta_markers, delta_markers, sex_gen_diet, file = "cell_type_markers.RData")



load("cell_type_markers.RData")

weight_16wk <- matched_phenos$weight_16wk
food_ave <- matched_phenos$food_ave
phenotypes <- cbind(sex_gen_diet, weight_16wk, food_ave)
head(phenotypes)

save(alpha_markers, beta_markers, delta_markers, phenotypes, file = "cell_type_markers2.RData")



























