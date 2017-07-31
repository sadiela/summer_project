# Genome Scans
# Sadie Allen
# July 31, 2017
# Genome scans of Alpha, Beta, Delta, Ins, BW16, and FoodAve w sex as add_cov, then int_cov

# Load libraries
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)

###
# prepare to run QTL scans

# load data 
load("/Users/s-allens/Documents/ssp/summer_project/data/DO378_islet.RData")
rownames(annot.samples) <- annot.samples$Mouse.ID

load("data/qtl_prep.RData")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"

int_covar <- model.matrix(~Sex, data = expr.adj)[,-1]
int_covar <- as.matrix(int_covar)
rownames(int_covar) <- annot.samples$Mouse.ID

# data frame rownames
expr.adj <- as.data.frame(expr.adj)
rownames(expr.adj) <- annot.samples$Mouse.ID


###
# Performing Scans
# Alpha, Beta, Delta, Ins, BW16 and FoodAvg
# addcovar
qtl.Alpha <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Alpha", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Beta <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Beta", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Delta <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Delta", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.Ins <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Ins", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.BW16 <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "BW16", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
qtl.FoodAvg <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "FoodAvg", drop = FALSE],
                   kinship = kin, addcovar = add_covar, cores = 4)
# intcovar
qtl.Alpha_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Alpha", drop = FALSE],
                   kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)
qtl.Beta_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Beta", drop = FALSE],
                  kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)
qtl.Delta_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Delta", drop = FALSE],
                   kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)
qtl.Ins_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "Ins", drop = FALSE],
                 kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)
qtl.BW16_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "BW16", drop = FALSE],
                  kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)
qtl.FoodAvg_int <- scan1(genoprobs = probs, pheno = expr.adj[,colnames(expr.adj) == "FoodAvg", drop = FALSE],
                     kinship = kin, addcovar = add_covar, intcovar = int_covar, cores = 4)

# differences
qtl.alpha_diff <- qtl.Alpha_int - qtl.Alpha
qtl.beta_diff <- qtl.Beta_int - qtl.Beta
qtl.delta_diff <- qtl.Delta_int - qtl.Delta
qtl.ins_diff <- qtl.Ins_int - qtl.Ins
qtl.bw16_diff <- qtl.BW16_int - qtl.BW16
qtl.foodavg_diff <- qtl.FoodAvg_int - qtl.FoodAvg

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.Alpha, map = map, main = "Alpha PC", col = "black")
plot_scan1(qtl.Alpha_int, map = map, main = "Alpha PC (int sex)", col = "black")
plot_scan1(qtl.alpha_diff, map = map, main = "Alpha diff", col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.Beta, map = map, main = "Beta PC", col = "black")
plot_scan1(qtl.Beta_int, map = map, main = "Beta PC (int sex)", col = "black")
plot_scan1(qtl.beta_diff, map = map, main = "Beta diff", col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.Delta, map = map, main = "Delta PC", col = "black")
plot_scan1(qtl.Delta_int, map = map, main = "Delta PC (int sex)", col = "black")
plot_scan1(qtl.delta_diff, map = map, main = "Delta diff", col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.Ins, map = map, main = "Ins PC", col = "black")
plot_scan1(qtl.Ins_int, map = map, main = "Ins PC (int sex)", col = "black")
plot_scan1(qtl.ins_diff, map = map, main = "Ins diff", col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.BW16, map = map, main = "BW16", col = "black")
plot_scan1(qtl.BW16_int, map = map, main = "BW16 (int sex)", col = "black")
plot_scan1(qtl.bw16_diff, map = map, main = "BW16 diff", col = "black")

quartz()
par(mfrow = c(3,1))
plot_scan1(qtl.FoodAvg, map = map, main = "FoodAvg", col = "black")
plot_scan1(qtl.FoodAvg_int, map = map, main = "FoodAvg (int sex)", col = "black")
plot_scan1(qtl.foodavg_diff, map = map, main = "FoodAvg diff", col = "black")






















