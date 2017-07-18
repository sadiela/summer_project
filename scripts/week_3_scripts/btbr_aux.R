# Week3
# Auxiliary data from BTBR data set: confirm findings of DiGruccio's paper?
# setwd("/Users/s-allens/Documents/ssp/summer_project")
library(qtl)

load("/Users/s-allens/Documents/ssp/summer_project/data/BTBR.clean.data.Rdata")

f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]

ghrl_exp <- islet.rz[, annot$a_gene_id[which(annot$gene_symbol=="Ghrl")]]
sst_exp <- islet.rz[, annot$a_gene_id[which(annot$gene_symbol=="Sst")]]
Sstr3_exp <- islet.rz[, annot$a_gene_id[which(annot$gene_symbol=="Sstr3")]]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum", "Sex", "pgm")], 
                   phenotypes.rz[c("Fat.wt", "GLU.10wk", "INS.10wk")], ghrl_exp, sst_exp, Sstr3_exp)

f2g<- calc.genoprob(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)
f2g <- sim.geno(f2g, step = 1, stepwidth = "fixed", map.function = "c-f", error.prob = 0.01)

sex <- as.numeric(f2g$pheno$Sex)

#f2g.perm1 <-scanone(f2g, pheno.col = 4:6, addcovar = sex, method= "hk", n.perm=100, perm.Xsp = TRUE)
#save(list="f2g.perm1", file="data/btbr_f2g_perm1.Rdata")

load(file="data/f2g_perm1.Rdata")

f2g.scanghrl <-scanone(f2g,  pheno.col = 7, addcovar = sex, method = "hk")

for(i in 1:6) {
  plot(f2g.scan1, lodcolumn = i)
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")
}

plot(f2g.scanghrl, lodcolumn = 1)
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

summary(f2g.scanghrl, chr=5)
#chr 5, pos 67.26

lodint(f2g.scanghrl, 5, 1.5, expandtomarkers = TRUE)
#position range for peak: 57.14493 - 79.39592 (rs13478458-rs13478542)


#Sema3a

sema3a_exp <- islet.rz[, annot$a_gene_id[which(annot$gene_symbol=="Sema3a")]]

f2g.scanghrl.sema3a <- scanone(f2g, pheno.col = 7, addcovar=c(sema3a_exp, Sex), method = "hk")
plot(f2g.scanghrl)
plot(f2g.scanghrl, f2g.scanghrl.sema3a, col=c("black", "red"))
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

f2g.scanins <-scanone(f2g,  pheno.col = 5, addcovar = sex, method = "hk")

f2g.scanins.sst <-scanone(f2g,  pheno.col = 5, addcovar =sst_exp, method = "hk")
plot(f2g.scanins, f2g.scanins.sst, lodcolumn = 1, col=c("black", "red"))
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")

f2g.scanins.sstr3 <- scanone(f2g,  pheno.col = 5, addcovar =sstr3_exp, method = "hk")
plot(f2g.scanins, f2g.scanins.sstr3, lodcolumn = 1, col=c("black", "red"))
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.05, lty="dashed", lwd=2, col="red")
add.threshold(f2g.scan1, perms=f2g.perm1, alpha=0.63, lty="dashed", lwd=2, col="green")


colnames(f2g$pheno)

















