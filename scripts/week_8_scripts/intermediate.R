# Find mediators of delta eigengene QTL peaks
# Sadie Allen
# July 26, 2017

QD18 <- get_genoprob(18, 5.990333)


source("Intermediate_Scripts/plot.mediation.R")
source("Intermediate_Scripts/mediation.scan.R")
source("Intermediate_Scripts/gmb.coordinates.R")
load("data/mouse.chrlen.rda")


med18 <- mediation.scan(target = delta_eigengene, mediator = rankz.mrna, annotation = annot.mrna, 
                      covar = add_covar, qtl.geno = QD18)
quartz()
plot(med18)  

med6 <- mediation.scan(target = delta_eigengene, mediator = rankz.mrna, annotation = annot.mrna, 
                      covar = add_covar, qtl.geno = QD6)
quartz()
plot(med)



anova(lm(delta_eigengene~sex + QD18))

anova(lm(delta_eigengene~sex + armc4_exp))

anova(lm(delta_eigengene~sex + armc4_exp + QD18))

anova(lm(delta_eigengene~sex + QD18 + armc4_exp))


anova(lm(delta_eigengene~sex + QD18))

anova(lm(delta_eigengene~sex + zfp438_exp))

anova(lm(delta_eigengene~sex + zfp438_exp + QD18))

anova(lm(delta_eigengene~sex + QD18 + zfp438_exp))



ins2_exp <- get_exp_dat("Ins2")
anova(lm(delta_eigengene~sex + ins2_exp))

anova(lm(delta_eigengene~sex + weight_16wk))


anova(lm(delta_eigengene~sex + ins2_exp + weight_16wk))
# effect is decreased

# Is body weight affecting alpha and delta eigengenes purely through increasing the number of beta cells?

anova(lm(delta_eigengene~sex + weight_16wk))
# p-value: 1.179e-12 ***
anova(lm(delta_eigengene~sex + ins2_exp + weight_16wk))
# does weight still have an effect on the eigengene after accounting for ins2_exp?
# YES!!!!
# p-value: 3.118e-05 ***, lower, but still a significant effect
# this tells us that body weight is still affecting delta cells in other ways

anova(lm(alpha_eigengene~sex + weight_16wk))
# p-value: < 2.2e-16 ***
anova(lm(alpha_eigengene~sex + ins2_exp + weight_16wk))
# p-value: 6.537e-08 ***

















