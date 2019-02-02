library ( emmeans )
library(lsmeans)
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
pigs.emm.s <- emmeans(pigs.lm, "source")

pairs(pigs.emm.s)
plot(pigs.emm.s, comparisons = TRUE, int.adjust = "none")



bd.lm = lm(value ~ X2 + X3 , data = bd)
bd.emm  <- emmeans(bd.lm, "X2")
pairs(bd.emm)
plot(bd.emm, comparisons = TRUE )






result = lsmeans::lsmeans(bd.lm,
        pairwise ~ X2,
        adjust="tukey")

r = summary ( result )
r$contrasts$p.value

c( "Control - MMD-hem" , "Control - MMD-ischemic", "MMD-hem - MMD-ischemic" )


bd$X2 = bd$X2 - bd$X3

bd.lm = lm(value ~ X2  , data = bd)
bd.emm  <- emmeans(bd.lm, "X2")
pairs(bd.emm)