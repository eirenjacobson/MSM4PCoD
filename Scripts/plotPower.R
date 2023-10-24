
id <- "D50AB_Real_2023-10-13"

plotPower <- function(id){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)


load(paste0("./Results/Trend_", id, ".RData"))


sim.results <- unique(trend.results) %>% filter(Type == "Sim") %>% select(-Sig)

sig.results <- unique(trend.results) %>% filter(Type != "Sim") 

tdf <- left_join(sim.results, sig.results, by = "Iteration", relationship = "one-to-many")

p1 <- ggplot(tdf) + 
  geom_smooth(data = filter(tdf, Type.y == "IPM"), aes(x=DeltaTrend.x, y = DeltaTrend.y), method = "lm", se = FALSE, lty = "dashed", color = "black") +
  geom_line(aes(x=DeltaTrend.x, y = DeltaTrend.x)) + # line = Sim
  geom_point(aes(x=DeltaTrend.x, y = DeltaTrend.y, color = factor(Sig), shape = factor(Type.y))) +
  xlab("Simulated Delta Trend")+
  ylab("Estimated Delta Trend")+
  theme_bw()

tdf <- unique(trend.results) %>% 
  pivot_wider(names_from = Type, values_from = c(DeltaTrend, Sig)) 

pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = filter(tdf, DeltaTrend_Sim<=0), family = "binomial")
pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = filter(tdf, DeltaTrend_Sim<=0), family = "binomial")
pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = filter(tdf, DeltaTrend_Sim<=0), family = "binomial")
newdata <- data.frame("DeltaTrend_Sim" = seq(-0.05, 0, by = 0.001))
newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")

nd <- newdata %>% pivot_longer(cols = c("IPM_Pred", "PAM_Pred", "LT_Pred"), names_to = "Type", values_to = "Pred")

p2 <- ggplot(nd) +
  geom_line(aes(x=abs(DeltaTrend_Sim), y = Pred, color = Type)) +
  theme_bw() +
  xlab("Simulated Annual Decline") +
  ylab("Predicted Power")+
  ylim(c(0,1))

ggsave(ggarrange(p1, p2, ncol = 2), 
       filename = paste0("./Figures/PowerPlots_", id, ".png"), 
       width = 12, height = 6, units = "in")

save(tdf, file = paste0("./Results/ObsPower_", id))
#save(nd, file = paste0("./Results/PredPower_", id))
}
