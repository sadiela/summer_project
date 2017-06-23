# Coat Color Example
# Sadie Allen
# June 23, 2017
# Using QTL2 to map the albino locus in mice

## Load data ##
# islet data
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID
# phenotype data
matched_phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)
ghrelin_shortlist <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/ghrelin_shortlist2.csv")
rownames(ghrelin_shortlist) <- ghrelin_shortlist[,1]
mod_ghrelin_shortlist <- cbind(mod_ghrelin_shortlist, coat_colors)

# Get coat color pheno
coat_colors <- matched_phenotypes$coat_color
mod_ghrelin_shortlist <- cbind(ghrelin_shortlist, coat_colors)


# Edit so there are only three colors: white, black, or agouti
coat_colors <- replace(coat_colors, coat_colors=="chinchilla", "agouti")
coat_colors <- replace(coat_colors, coat_colors == "WSB", "agouti")

# Numerically assign coat colors
coat_colors <- replace(coat_colors, coat_colors=="white", "1")
coat_colors <- replace(coat_colors, coat_colors=="black", "-1")
coat_colors <- replace(coat_colors, coat_colors=="agouti", "0")
coat_colors <- as.numeric(coat_colors)

# Convert genoprobs to QTL format
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")

# This just changes the chromosome names "X" to "20"
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

# Calculate kinship
kin <- calc_kinship(probs = probs, type = "loco", cores = 4)

# Create additive covariate
temp = merge(annot.samples, ghrelin_shortlist, by = "row.names")
rownames(temp) = temp[,1]
# annot.samples <- merge(annot.samples, diet_days_id, by.x="Mouse.ID", by.y = "Mouse.ID")
# Now, create the covariate
add_covar <- model.matrix(~Sex + Generation + diet_days, data = temp)[,-1]

#convert a marker map organized as data frame to a list
map <- map_df_to_list(map = snps, pos_column = "bp")

# Do scans
qtl.coat_color <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,24, drop = FALSE],
                        cores = 4)

qtl.coat_color.wkin <- scan1(genoprobs = probs, pheno = mod_ghrelin_shortlist[,24, drop = FALSE],
                        kinship = kin, cores = 4)

plot_scan1(x = qtl.coat_color, map = map, main = colnames(qtl.coat_color)[1])
plot_scan1(x = qtl.coat_color.wkin, map = map, main = colnames(qtl.coat_color.wkin)[1])

chr = 7
qtl.coat_color.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,24, drop = FALSE])
plot(x = qtl.coat_color.coef, scan1_output = qtl.coat_color, map = map[[chr]], columns = 1:8, col = CCcolors)

chr = 7
qtl.coat_color.wkin.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,24, drop = FALSE], kinship = kin[[chr]])
plot(x = qtl.coat_color.coef, scan1_output = qtl.coat_color.wkin, map = map[[chr]], columns = 1:8, col = CCcolors)

chr = 2
qtl.coat_color.coef <- scan1coef(genoprobs = probs[,chr], pheno = mod_ghrelin_shortlist[,24, drop = FALSE])
plot(x = qtl.coat_color.coef, scan1_output = qtl.coat_color, map = map[[chr]], columns = 1:8, col = CCcolors)


# This seemed to work. The models with and without kinship were almost, if not completely, iidentical,
# which suggests that the kinship matrix is not causing problems with my other scans. This leads
# me to believe I might have done something to change my phenotypes. 



