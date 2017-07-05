# Simple BIC Modeling
# Sadie Allen

cor(ghsr_exp, ghrelin_list$food_ave)
cor(svil_exp, ghsr_exp)
Q18 <- get_genoprob(18, 4.474313)

efnb3_exp <- get_exp_dat("Efnb3")
cor(ghrelin_list$food_ave, efnb3_exp)
ghrelin_list <- cbind(ghrelin_list, ghsr_exp, efnb3_exp, Q18)

ggplot(ghrelin_list, aes(x = ghsr_exp, y = efnb3_exp, group = sex, fill = sex, col = sex)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm", formula = y~x)

ggplot(ghrelin_list, aes(x = food_ave, y = efnb3_exp, group = sex, fill = sex, col = sex)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm", formula = y~x)

ggplot(ghrelin_list, aes(x = svil_exp, y = ghsr_exp, group = sex, fill = sex, col = sex)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm", formula = y~x)

ggplot(ghrelin_list, aes(x = svil_exp, y = food_ave, group = sex, fill = sex, col = sex)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm", formula = y~x)

# Sue really wants me to try simple BIC modeling so I will test the Svil relationship
# and the Efnb3 relationship

# Svil expression affects ghsr expression
BIC(lm(ghsr_exp~ghrelin_list$sex))
# 1009.857
BIC(lm(ghsr_exp~ghrelin_list$sex + svil_exp))
#997.6577

#Q18 affects svil expression
BIC(lm(svil_exp~ghrelin_list$sex))
#1065.1
BIC(lm(svil_exp~ghrelin_list$sex + Q18))
#916.3928

#Q18 affects ghsr expression
BIC(lm(ghsr_exp~ghrelin_list$sex))
# 1009.857
BIC(lm(ghsr_exp~ghrelin_list$sex + Q18))
#1009.763
BIC(lm(ghsr_exp~ghrelin_list$sex + svil_exp + Q18))
#1011.416




BIC(lm(ghrelin_list$food_ave~ghrelin_list$sex))
# 472.1798
BIC(lm(ghrelin_list$food_ave~ghrelin_list$sex + ghsr_exp))
# 423.0559
BIC(lm(ghrelin_list$food_ave~ghrelin_list$sex + ghsr_exp + efnb3_exp))
# 409.1854













