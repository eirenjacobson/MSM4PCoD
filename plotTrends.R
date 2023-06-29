
library(dplyr)
library(tidyr)
library(ggplot2)

load("./Results/Trend_D50_LCP_2023-06-27.RData")

tdf <- trend.results %>% 
  pivot_wider(names_from = Type, values_from = c(Trend, Sig)) %>%
  filter(Period == "Post")

ggplot(tdf) + 
  geom_point(aes(x=Trend_Sim, y = Trend_Sim, fill = factor(Sig_Sim)), pch = 21) + # circle = Sim
  geom_point(aes(x= Trend_Sim, y = Trend_LT,fill = factor(Sig_LT,)), pch = 22) + # square = LT
  geom_point(aes(x=Trend_Sim, y = Trend_PAM, fill = factor(Sig_PAM)), pch = 23)+ # diamond = PAM
  geom_point(aes(x= Trend_Sim, y = Trend_Model, fill = factor(Sig_IPM)), pch = 24)+ # triangle = IPM
  xlab("Simulated Trend")+
  ylab("Estimated Trend")+
  theme_bw()