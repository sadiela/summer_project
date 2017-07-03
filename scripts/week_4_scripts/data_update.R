# Update to New Phenotype Data
# Sadie Allen
# June 29, 2017
# Updating to new completed phenotype data provided by Karl Broman

# Load files from the email
load("/Users/s-allens/Library/Caches/TemporaryItems/Outlook Temp/pheno_clin.RData")

# Save the phenotypes as RData
save(pheno_clin, file="/Users/s-allens/Documents/ssp/summer_project/data/pheno_clin_update.RData")

# Save dictionary as a csv file
write.csv(pheno_clin_dict, file="/Users/s-allens/Documents/ssp/summer_project/data/pheno_dict_update")

# Load islet data for matching
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")

# Now must match phenotyped mice with genotyped mice
# Load neccessary library:
library(dplyr)

# Load functions
source("scripts/functions.R")

# add 0 to the identifiers for mouse ids < 100
pheno_clin$mouse[nchar(pheno_clin$mouse) == 5] <- sub("DO-", "DO-0", pheno_clin$mouse[nchar(pheno_clin$mouse) == 5])

#there are more mice in the phenotype data than in the mRNA data, we need to match up the ones
#that are in both and delete the extras
names(pheno_clin)[1] <- "Mouse.ID"
idx <- match(annot.samples$Mouse.ID, pheno_clin$Mouse.ID)

# Cuts down pheno dataset to 378 animals, the same ones for which RNAseq data is provided
pheno_clin <- pheno_clin[idx,]

# Check that the IDs for mRNA data and phenotype data are the same
stopifnot(annot.samples$Mouse.ID == pheno_clin$Mouse.ID)

# Change variable names
matched_phenos <- pheno_clin

# Set row names
rownames(matched_phenos) <- matched_phenos$Mouse.ID

# Save new matched phenotype data
write.csv(matched_phenos, file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv")
# TO READ IN:
matched_phenos <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_phenos.csv", as.is=TRUE)

# Now, get list of relevant phenotypes for ghrelin project
ghrelin_list <- matched_phenos %>% select(Mouse.ID, sex, diet_days, Glu_0min,  Glu_sac, Ins_0min, 
                         Ins_sac, food_ave, weight_sac)

# Get GR expression data and put it into the list
ghsr_exp <- get_exp_dat("Ghsr")
ghrelin_list <- merge(ghrelin_list, ghsr_exp, by="row.names")

# Clean up list
names(ghrelin_list)[11] <- "ghsr_exp"
ghrelin_list$Row.names <- NULL

# Save this list
write.csv(ghrelin_list, file="/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")
ghrelin_list <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_list.csv")

