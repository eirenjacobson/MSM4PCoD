
# compare power with and without calf obs submodel

library(dplyr)
library(tidyr)
library(ggplot2)

# ids: vector of scenarios to be combined
# type: added to filename as indicator

#type <- "Opt" 

ids_T1 <- c("NULL_Real_2023-10-13", 
         "D50AB_Real_2023-10-13")

ids_T2 <- c("NULL_Real_wCalfData_2023-10-10",
            "D50AB_Real_wCalfData_2023-10-12")

nfiles <- length(ids_T1)

data_T1 <- data.frame()

for (i in 1:nfiles){
  
  load(paste0("./Results/ObsPower_", ids_T1[i]))
  tdf$id <- ids_T1[i]
  data_T1 <- rbind(data_T1, tdf)
  
} # end for i


pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = filter(data_T1, DeltaTrend_Sim<=0), family = "binomial")
pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = filter(data_T1, DeltaTrend_Sim<=0), family = "binomial")
pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = filter(data_T1, DeltaTrend_Sim<=0), family = "binomial")
newdata <- data.frame("DeltaTrend_Sim" = seq(range(data_T1$DeltaTrend_Sim)[1], 0, by = 0.001))
newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")

nd_1 <- newdata %>% 
  pivot_longer(cols = c("IPM_Pred", "PAM_Pred", "LT_Pred"), names_to = "Type", values_to = "Pred") %>%
  mutate(CalfRatio = FALSE)


nfiles <- length(ids_T2)

data_T2 <- data.frame()

for (i in 1:nfiles){
  
  load(paste0("./Results/ObsPower_", ids_T2[i]))
  tdf$id <- ids_T2[i]
  data_T2 <- rbind(data_T2, tdf)
  
} # end for i


pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = filter(data_T2, DeltaTrend_Sim<=0), family = "binomial")
pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = filter(data_T2, DeltaTrend_Sim<=0), family = "binomial")
pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = filter(data_T2, DeltaTrend_Sim<=0), family = "binomial")
newdata <- data.frame("DeltaTrend_Sim" = seq(range(data_T2$DeltaTrend_Sim)[1], 0, by = 0.001))
newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")

nd_2 <- newdata %>% 
  pivot_longer(cols = c("IPM_Pred", "PAM_Pred", "LT_Pred"), names_to = "Type", values_to = "Pred") %>%
  mutate(CalfRatio = TRUE)





nd <- rbind.data.frame(nd_1, nd_2)



p2 <- ggplot(filter(nd, Type == "IPM_Pred")) +
  geom_line(aes(x=abs(DeltaTrend_Sim), y = Pred, color = CalfRatio)) +
  theme_bw() +
  xlab("Simulated Annual Decline") +
  ylab("Predicted Power")+
  ylim(c(0,1))

ggsave(p2, 
       filename = paste0("./Figures/ComparePowerPlots_wCalfRatio", type, ".png"), 
       width = 7, height = 6, units = "in")
