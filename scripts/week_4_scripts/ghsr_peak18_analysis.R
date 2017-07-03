# Ghsr Peak on Chromosome 18
# Sadie Allen
# June 27, 2017
# QTL scan for Ghsr


## Step 1: Find QTL peak for ghsr (14) ##
# GARY QUESTION: What is the threshold for a significant QTL peak? 6, 8, 11

###### QTL PREP ######## (needed for all qtl scans), don't have to run this anymore, just load RData file

load("data/qtl_prep.RData")

# Convert genoprobs from DOQTL to QTL2 format
#probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Calculate kinship
#kin <- calc_kinship(probs = probs, type = "loco", cores = 4)

# Additive covariates
#temp = merge(annot.samples, ghrelin_shortlist, by = "row.names")
#rownames(temp) = temp[,1]
# annot.samples <- merge(annot.samples, diet_days_id, by.x="Mouse.ID", by.y = "Mouse.ID")
# Now, create the covariate
#add_covar <- model.matrix(~Sex + Generation + diet_days, data = temp)[,-1]

# Convert a marker map organized as data frame to a list
#map <- map_df_to_list(map = snps, pos_column = "bp")

### RUNNING AND PLOTTING SCANS ###

# Perform scan
qtl.ghsr <- scan1(genoprobs = probs, pheno = ghrelin_list[,colnames(ghrelin_list)=="ghsr_exp", drop = FALSE],
                  kinship = kin, addcovar = add_covar, cores = 4)
# Plot scan
plot_scan1(x = qtl.ghsr, map = map, main = colnames(qtl.ghsr)[1])
# Find peak positions
find_peaks(qtl.ghsr, map, threshold = 6, drop = 1.5) # What does the drop argument do?
# chromosome 4, pos 11.984686, lod = 6.859136 (interval 0.207518 - 13.32408)
# chromosome 18, pos 4.474313, lod = 7.822138 (interval 0.026149 - 10.21129)

Nhs_exp <- get_exp_dat("Nhs")
ghrelin_shortlist <- cbind(ghrelin_shortlist, Nhs_exp)
qtl.Nhs <- scan1(genoprobs = probs, pheno = ghrelin_shortlist[,colnames(ghrelin_shortlist)=="Nhs_exp", drop = FALSE],
                 kinship = kin, addcovar = add_covar, cores = 4)
plot_scan1(x = qtl.Nhs, map = map, main = colnames(qtl.Nhs)[1])
find_peaks(qtl.Nhs, map, threshold = 6, drop = 1.5) # What does the drop argument do?
# chromosome 18, pos 43.10387, LOD 9.471 (located on the X chromosome...)
###################################################################################################
