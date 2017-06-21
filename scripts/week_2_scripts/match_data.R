#Week2
#this script gets rid of the phenotypical data for mice that do not have gene
#expression data
#load phenotypical data from large .csv file
phenotypes <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/pheno_clin.csv", as.is=TRUE)

# for mouse id < 100, add 0 to the identifiers
phenotypes$mouse[nchar(phenotypes$mouse) == 5] <- sub("DO-", "DO-0", phenotypes$mouse[nchar(phenotypes$mouse) == 5])

#there are more mice in the phenotype data than in the mRNA data, we need to match up the ones
#that are in both and delete the extras
names(phenotypes)[1] <- "Mouse.ID"
idx <- match(annot.samples$Mouse.ID, phenotypes$Mouse.ID)

#cuts down phenotype dataset to 378 animals, the same ones for which RNAseq data is provided
phenotypes <- phenotypes[idx,]

#check that the IDs for mRNA data and phenotype data are the same
stopifnot(annot.samples$Mouse.ID == phenotypes$Mouse.ID)

#changing variable name
matched_phenotypes <- phenotypes

#This section works
#saving matched phenotypes to a new data file
write.csv(matched_phenotypes, file = "/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", row.names = FALSE)
#NOW THAT THIS HAS BEEN DONE ONCE, SIMPLY RUN THE FOLLOWING TO GET ACCESS TO THE MATCHED PHENOTYPES
matched_phenotypes2 <- read.csv("/Users/s-allens/Documents/ssp/summer_project/data/matched_pheno_clin.csv", as.is=TRUE)
