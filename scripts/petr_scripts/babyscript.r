#mouseIDs from mRNAseq data (annot.samples$Mouse.ID)
samplenums_mRNA <- annot.samples$Mouse.ID
#mouseIDs from phenotype data
samplenums_mRNA <- phenotypes$mouse

#create empty vectors to fill with the id numbers
mrna_ids <- character(0)
pheno_ids <- character(0)

for(i in phenotypes$mouse) {
  print(i)
  x <- substr(i, 4, nchar(i))
  print(x)
  pheno_ids <- c(pheno_ids, x)
}
pheno_ids

for(i in annot.samples$Mouse.ID) {
  #print(i)
  x <- substr(i, 4, nchar(i))
  #print(x)
  mrna_ids <- c(mrna_ids, x)
}
mrna_ids

