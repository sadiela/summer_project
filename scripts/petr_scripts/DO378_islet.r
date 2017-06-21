load("DO381_islet_emase_m4.RData")

env <- new.env()
load("DO_islet_gigamuga_grid.RData", envir=env)

# rename mice
pheno$mouse <- pheno[,6]
annot.samples <- pheno[,1:5] 
names(annot.samples) <- c("Mouse.ID", "Sex", "Generation", "Age", "Sample.Number")
annot.samples$Age <- as.numeric(annot.samples$Age)
# we filled G20 because EMASE crashes on G21, now it is time to correct for it
annot.samples$Generation[annot.samples$Generation=="G20"] <- "G21"
rm(pheno)

# add information about chrM and chrY
chrMY <- read.csv("JAX_gigamuga_master_ChrMY_v3.csv", as.is=TRUE)
chrMY <- chrMY[grep("Attie", chrMY$file),]
chrMY$sample.number <- as.numeric(sub(".*DO([0-9]*).*", "\\1", chrMY$sample.id))
idx <- match(annot.samples$Sample.Number, chrMY$sample.number)
annot.samples$chrM <- chrMY$chrM[idx]
annot.samples$chrY <- chrMY$chrY[idx]


# reorder GigaMUGA data
muga.number <- as.numeric(sub(".*DO([0-9]*).*", "\\1", dimnames(env$probs)[[1]]))
idx <- match(annot.samples$Sample.Number, muga.number)
mugaprobs <- env$probs[idx,,]

cors <- rep(NA, nrow(probs)) 
for (i in 1:nrow(probs)) {
  print(i)
  cors[i] <- cor(as.vector(probs[i,,]), as.vector(mugaprobs[i,,])) 
}
summary(cors)

# suspicious samples
annot.samples[which(cors < 0.6),]
#                 Mouse.ID Sex Generation Age Sample.Number
# 248-GES15-05892   DO-248   F        G19  19           248
# 268-GES15-05657   DO-268   F        G19  21           268
# 269-GES15-05831   DO-269   F        G19  21           269
# 310-GES15-05719   DO-310   M        G19  20           310
# 360_GES16_02094   DO-360   F        G21  19           360
# 370_GES16_02127   DO-370   F        G21  20           370

# for samples that are not suspicious, copy haplotype probs
for (i in which(cors >= 0.6)) {
  probs[i,,] <- mugaprobs[i,,]
}

library(HPQTL)
p2 <- probs
attr(p2, "markers") <- snps
geno <- extract.geno(p2)
geno$subjects <- annot.samples[,1]; 
geno$calls <- LETTERS[1:8]
Glist <-  gensim.matrix(geno, procedure="LOCO", verbose=TRUE)
for (i in 1:length(Glist))
  rownames(Glist[[i]]) <- colnames(Glist[[i]]) <- annot.samples[,1]
G <- gensim.matrix(geno, procedure="LMM", verbose=TRUE)

genoprobs <- probs
rm(probs)

covar <- model.matrix(~Sex+Generation, data=annot.samples)

colnames(G) <- rownames(G) <- annot.samples[,1]
for (i in 1:20) colnames(Glist[[i]]) <- rownames(Glist[[i]]) <- annot.samples[,1]
rownames(expr.mrna) <- rownames(raw.mrna) <- rownames(covar) <- annot.samples[,1]
dimnames(genoprobs)[[1]] <- annot.samples[,1]

# remove duplicated samples and outlier DO-248 (by Aimee Broman)
to.be.removed <- which(duplicated(annot.samples$Mouse.ID) | annot.samples$Mouse.ID == "DO-248")
stopifnot(length(to.be.removed) == 3)
genoprobs <- genoprobs[-to.be.removed,,]
expr.mrna <- expr.mrna[-to.be.removed,]
raw.mrna <- raw.mrna[-to.be.removed,]
G <- G[-to.be.removed,-to.be.removed]
for (i in 1:length(Glist))
  Glist[[i]] <- Glist[[i]][-to.be.removed,-to.be.removed]
annot.samples <- annot.samples[-to.be.removed,]
covar <- covar[-to.be.removed,] 

## combat normalization
library(Biobase)
library(sva)

# create BioC ExpressionSet object
texpr.mrna <- t(expr.mrna)
colnames(texpr.mrna) <- make.names(colnames(texpr.mrna), unique=TRUE)
row.names(annot.samples) <-  colnames(texpr.mrna)
row.names(annot.mrna) <- row.names(texpr.mrna)
phenodata <- new("AnnotatedDataFrame", data=annot.samples)
featuredata <- new("AnnotatedDataFrame", data=annot.mrna)
eset <- ExpressionSet(assayData=texpr.mrna, phenoData=phenodata, featureData=featuredata)
edata = exprs(eset)

# apply ComBat and save the results
mod0 = model.matrix(~Sex,data=annot.samples)
combat_edata = ComBat(dat=edata, batch=annot.samples$Generation, mod=mod0, par.prior=TRUE, prior.plots=FALSE)
combat.mrna <- t(combat_edata)
rownames(combat.mrna) <- rownames(expr.mrna)

# switch names
rankz.mrna <- expr.mrna # expression after rankZ normalization
expr.mrna <- combat.mrna # after rankZ and ComBat

# for Matt's browser
covar_factors <- data.frame(column_name=c("Sex", "Generation"), 
         display_name=c("Sex","Generation"), stringsAsFactors=FALSE)

# fix for row-names of annot.samples
rownames(annot.samples) <- rownames(covar) <- rownames(expr.mrna)

save(genoprobs, annot.mrna, expr.mrna, rankz.mrna, raw.mrna, covar, G, Glist, annot.samples, snps, covar_factors, file="DO378_islet.RData")

# cisQTL check
lods <- function(expr, prob, id, X) {
  results <- rep(0,ncol(expr))
  for(i in which(!is.na(id))) {
    if (i %% 1000 ==0) print(i)
    y <- expr[,i] # expression of one gene/protein
    sel <- !is.na(y)
    rss1 <- sum(lsfit(y=y[sel], x=cbind(X[sel,],prob[sel,-1,id[i]]), intercept=FALSE)$residuals^2)
    rss0 <- sum(lsfit(y=y[sel], x=X[sel,], intercept=FALSE)$residuals^2)
    results[i] <- sum(sel)/2 * (log10(rss0) - log10(rss1))
  }
  results
}

qgrid <- seq(from=0.75, to=0.95, length=100)
tmp <- lods(rankz.mrna, genoprobs, annot.mrna$nearest_snp, covar)
print(quantile(tmp, qgrid, na.rm=TRUE))

tmp <- lods(expr.mrna, genoprobs, annot.mrna$nearest_snp, covar)
print(quantile(tmp, qgrid, na.rm=TRUE))


