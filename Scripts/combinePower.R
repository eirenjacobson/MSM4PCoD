# Script to plot power across a broader range of DeltaTrend
# by combining multiple scenarios (e.g. NULL + a decline)

combinePower <- function(ids, type){
  
  # ids: vector of scenarios to be combined
  # type: added to filename as indicator
  
  type <- "Real" 
  
  ids <- c("NULL_Real_2023-10-13", 
           "D50AB_Real_2023-10-13")
  
  nfiles <- length(ids)
  
  data <- data.frame()
  
  for (i in 1:nfiles){
    
    load(paste0("./Results/ObsPower_", ids[i]))
    tdf$id <- ids[i]
    data <- rbind(data, tdf)
    
  } # end for i
  
  
  pow_ipm <- glm(Sig_IPM ~ DeltaTrend_Sim, data = filter(data, DeltaTrend_Sim<=0), family = "binomial")
  pow_pam <- glm(Sig_PAM ~ DeltaTrend_Sim, data = filter(data, DeltaTrend_Sim<=0), family = "binomial")
  pow_lt <- glm(Sig_LT ~ DeltaTrend_Sim, data = filter(data, DeltaTrend_Sim<=0), family = "binomial")
  newdata <- data.frame("DeltaTrend_Sim" = seq(range(data$DeltaTrend_Sim)[1], 0, by = 0.001))
  newdata$IPM_Pred <- predict(pow_ipm, newdata = newdata, type = "response")
  newdata$PAM_Pred <- predict(pow_pam, newdata=newdata, type = "response")
  newdata$LT_Pred <- predict(pow_lt, newdata=newdata, type = "response")
  
  nd <- newdata %>% pivot_longer(cols = c("IPM_Pred", "PAM_Pred", "LT_Pred"), names_to = "Type", values_to = "Pred")
  
  p2 <- ggplot(nd) +
    geom_line(aes(x=DeltaTrend_Sim, y = Pred, color = Type)) +
    theme_bw() +
    ylim(c(0,1))
  
  ggsave(p2, 
         filename = paste0("./Figures/CombinedPowerPlots_", type, ".png"), 
         width = 7, height = 6, units = "in")
  
  
  tdata <- data.frame("DeltaTrend_Sim" = c(0, -0.005, -0.01, -0.02))
  tdata$IPM_Pred <- predict(pow_ipm, newdata = tdata, type = "response")
  tdata$PAM_Pred <- predict(pow_pam, newdata=tdata, type = "response")
  tdata$LT_Pred <- predict(pow_lt, newdata=tdata, type = "response")
  
  round(tdata, digits = 2)
} # end function

