

load("./Results/ObsPower_D50B_Real_wCalfData_2023-08-21")
mega.tdf <- tdf
load("./Results/ObsPower_D75B_Real_wCalfData_2023-08-21")
mega.tdf <- rbind.data.frame(tdf, mega.tdf)
load("./Results/ObsPower_NULL_Real_wCalfData_2023-08-21")
mega.tdf <- rbind.data.frame(mega.tdf, tdf)



pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = mega.tdf, family = "binomial")
pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = mega.tdf, family = "binomial")
pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = mega.tdf, family = "binomial")
newdata <- data.frame("DeltaTrend_Sim" = seq(-0.08, 0.01, by = 0.001))
newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")

nd <- newdata %>% pivot_longer(cols = c("IPM_Pred", "PAM_Pred", "LT_Pred"), names_to = "Type", values_to = "Pred")

p2 <- ggplot(nd) +
  geom_line(aes(x=DeltaTrend_Sim, y = Pred, color = Type)) +
  theme_bw() +
  ylim(c(0,1))

ggsave(plot = p2, filename = paste0("./Figures/PowerPlots_Real_wCalfData.png"), height = 5, width = 6)

