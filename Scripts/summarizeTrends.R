
library(tidyr)
library(dplyr)
library(ggplot2)

id <- 

summarizeTrends <- function(){
  
  
  load("./Results/Trend_D75_LCP_MidCval_2023-07-20.RData")
  
  #TODO plot this as simulated trend x axis, 
  tdf <- unique(trend.results) %>% 
    #select(-Sig) %>%
  #  group_by(c(Iteration, Type)) %>%
    #pivot_wider(names_from = Period, values_from = Trend) %>%
    #mutate("Delta" = Post - Pre) %>%
    filter(Type != "Sim") %>%
    filter(Type != "Model")
  

  ggplot(tdf) + 
    geom_density(aes(x=DeltaTrend))+
    facet_wrap(~Type, scales = "free_y") +
    theme_bw()
  
  ggplot(tdf) + 
    geom_density(aes(x=DeltaTrend))+
    facet_wrap(~Type, scales = "free") +
    theme_bw()
  
  
  
  widedf <- tdf %>% 
    group_by(Iteration) %>%
    pivot_wider(names_from = Type, values_from = c(DeltaTrend)) 
  
  
  ggplot(widedf) +
    geom_point(aes(x=Sim, y = Sim), pch = 21) + # circle = Sim
    geom_point(aes(x= Sim, y = LT), pch = 22) + # square = LT
    geom_point(aes(x=Sim, y = PAM), pch = 23)+ # diamond = PAM
    geom_point(aes(x= Sim, y = IPM), pch = 24)+ # triangle = IPM
    xlab("Simulated Trend")+
    ylab("Estimated Trend")+
    theme_bw()
  
  ggplot(widedf) +
    geom_point(aes(x=Sim, y = Sim), pch = 21) + # circle = Sim
    geom_point(aes(x= Sim, y = LT), pch = 22) + # square = LT
    geom_point(aes(x=Sim, y = PAM), pch = 23)+ # diamond = PAM
    geom_point(aes(x= Sim, y = IPM), pch = 24)+ # triangle = IPM
    xlab("Simulated Trend")+
    ylab("Estimated Trend")+
    theme_bw()
    
}