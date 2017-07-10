# Full Genome Mediation Analysis
# Sadie Allen
# July 10, 2017
# In this script I will run the expression of all 20,000+ genes through to check if they mediate 
# several different relationships.

rm(list =ls())

# Load data
# Phenotype Data
ghrelin_list <- read.csv(file = "/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
rownames(ghrelin_list) <- ghrelin_list$Mouse.ID
# Islet Data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

# deleting a row that has an obsolete gene symbol (there are two, deleting the old one)
# causing errors in mediation
annot.mrna <- annot.mrna[-16503, ]
annot.mrna <- annot.mrna[-9377,]
annot.mrna <- annot.mrna[-13944, ]

# Load functions
source("scripts/functions.R")

gene_names <- get_gene_names()

weight_sac <- ghrelin_list$weight_sac
food_ave <- ghrelin_list$food_ave
ghsr_exp <- ghrelin_list$ghsr_exp
sex <- ghrelin_list$sex

# First: look for mediators of the relationship between food intake and body weight
food_weight_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  gene <- gene_names[i]
  print(i)
  if(length(get_exp_dat(gene)) == 378) {
    if(mediation_analysis(food_ave, gene_names[i], weight_sac) == TRUE) {
      food_weight_mediators <- c(food_weight_mediators, gene_names[i])
    }
  }
}
# NO MEDIATORS!


# Next: look for gene that affects body weight and mediated through food intake
pre_food_weight_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  if(length(get_exp_dat(gene_names[i])) == 378) {
    if(first_mediation_analysis(gene_names[i], food_ave, weight_sac) == TRUE) {
      pre_food_weight_mediators <- c(pre_food_weight_mediators, gene_names[i])
    }
  }
}
pre_food_weight_mediators

# Look for mediators between ghsr expression and food_ave
ghsr_food_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  print(i)
  if(length(get_exp_dat(gene_names[i])) == 378 & gene_names[i] != "Ghsr") {
    if(mediation_analysis(ghsr_exp, gene_names[i], food_ave) == TRUE) {
      ghsr_food_mediators <- c(ghsr_food_mediators, gene_names[i])
    }
  }
}
ghsr_food_mediators
# Efnb3

# Next: look for gene that affects food intake and is mediated through ghsr expression
pre_ghsr_food_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  print(i)
  if(length(get_exp_dat(gene_names[i])) == 378 & gene_names[i] != "Ghsr") {
    if(first_mediation_analysis(gene_names[i], ghsr_exp, food_ave) == TRUE) {
      pre_ghsr_food_mediators <- c(pre_ghsr_food_mediators, gene_names[i])
    }
  }
}
pre_ghsr_food_mediators

# NOW WITH WEIGHT CHANGE
weight_change <- ghrelin_list$weight_change
food_weightchange_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  gene <- gene_names[i]
  print(i)
  if(length(get_exp_dat(gene)) == 378) {
    if(mediation_analysis(food_ave, gene_names[i], weight_change) == TRUE) {
      food_weightchange_mediators <- c(food_weightchange_mediators, gene_names[i])
    }
  }
}
# NO MEDIATORS


# Next: look for gene that affects body weight and mediated through food intake
pre_food_weightchange_mediators <- character()
for(i in 1:length(gene_names)) {
  print(gene_names[i])
  if(length(get_exp_dat(gene_names[i])) == 378) {
    if(first_mediation_analysis(gene_names[i], food_ave, weight_change) == TRUE) {
      pre_food_weightchange_mediators <- c(pre_food_weightchange_mediators, gene_names[i])
    }
  }
}
pre_food_weightchange_mediators

save(pre_food_weight_mediators, ghsr_food_mediators, pre_ghsr_food_mediators, pre_food_weightchange_mediators, file = "data/relationship_mediators.RData")





