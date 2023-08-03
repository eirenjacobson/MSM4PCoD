
library(dplyr)
library(tidyr)
library(ggplot2)

load("./Results/Trend_D50_LCP_Ideal1_2023-07-27.RData")

tdf <- unique(trend.results) %>% 
  pivot_wider(names_from = Type, values_from = c(DeltaTrend, Sig)) 

ggplot(tdf) + 
  geom_point(aes(x=DeltaTrend_Sim, y = DeltaTrend_Sim, fill = factor(Sig_Sim)), pch = 21) + # circle = Sim
#  geom_point(aes(x= DeltaTrend_Sim, y = DeltaTrend_LT,fill = factor(Sig_LT,)), pch = 22) + # square = LT
#  geom_point(aes(x=DeltaTrend_Sim, y = DeltaTrend_PAM, fill = factor(Sig_PAM)), pch = 23)+ # diamond = PAM
  geom_point(aes(x= DeltaTrend_Sim, y = DeltaTrend_IPM, fill = factor(Sig_IPM)), pch = 24)+ # triangle = IPM
  xlab("Simulated Delta Trend")+
  ylab("Estimated Delta Trend")+
  theme_bw()


pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = tdf, family = "binomial")
pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = tdf, family = "binomial")
pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = tdf, family = "binomial")
newdata <- data.frame("DeltaTrend_Sim" = seq(-0.05, 0, by = 0.001))
newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")

plot(newdata$DeltaTrend_Sim, newdata$PAM_Pred, type = "l", col = "blue", 
     ylim = c(0,1), ylab = "P(Change in Trend Detected)", xlab = "Simulated Change in Trend")
lines(newdata$DeltaTrend_Sim, newdata$LT_Pred, type = "l", col = "red")
lines(newdata$DeltaTrend_Sim, newdata$IPM_Pred, type = "l", col = "green")
