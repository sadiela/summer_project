# BTBR Ghrs connections
# Sadie Allen
# July 17, 2017
# Checking if Ghsr expression in the islet cells is related to its expression
# in the brain

rm(list = ls())

load("/Users/s-allens/Documents/ssp/summer_project/data/BTBR.clean.data.Rdata")


f2g$pheno
colnames(phenotypes.rz)

ghsr.islet <- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ghsr")]]
efnb3.islet <- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Efnb3")]]

ghsr.hypo <- hypo.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ghsr")]]
efnb3.hypo <- hypo.rz[,annot$a_gene_id[which(annot$gene_symbol=="Efnb3")]]

ghrl.islet <- islet.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ghrl")]]
ghrl.gastroc <- gastroc.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ghrl")]]
ghrl.hypo <- hypo.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ghrl")]]


f2g$pheno <- cbind(f2g$pheno[,c("MouseNum", "Sex", "pgm")],phenotypes.rz[,c("Weight", "WT.6wk")], 
                   ghsr.islet, ghsr.hypo, efnb3.islet, efnb3.hypo, ghrl.islet, ghrl.gastroc, ghrl.hypo)

names(f2g$pheno)

pcor <- cor(f2g$pheno[,4:12], use= "complete.obs")
round(pcor, digits=2)
corrplot(pcor )









