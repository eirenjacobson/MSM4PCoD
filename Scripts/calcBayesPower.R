
calcBayesPower <- function(id, pars){
  
  library(dplyr)
  library(tidyr)

  id <- "D50AB_Ideal_wCalfData_2023-10-06"
  #nsim <- 100
  #nsim <- pars$nsim
  
  load(paste0("./Results/ProcResults_", id, ".RData")) #results.out
  load(paste0("./Results/Trend_", id, ".RData")) #results.out

  
  intervals <- results.out$pardf %>% 
    group_by(Iter) %>%
    summarize(K1_LCI = quantile(K1, 0.025), 
              K2_LCI = quantile(K2, 0.025),
              K1_UCI = quantile(K1, 0.975),
              K2_UCI = quantile(K2, 0.975)) %>%
    rename(Iteration = Iter)
  
  intervals$BSig <- ifelse(intervals$K1_UCI < intervals$K2_LCI | intervals$K2_UCI < intervals$K1_LCI, TRUE, FALSE)
  
  sim.results <- unique(trend.results) %>% filter(Type == "Sim") %>% select(-Sig)
  
  tr <- filter(trend.results, Type == "IPM") %>% left_join(intervals, by = "Iteration")
               
  neg <- filter(tr, DeltaTrend<0)
  pos <- filter(tr, DeltaTrend>0)
  pow1n_model <- glm(Sig ~ DeltaTrend, data = neg, family = "binomial")
  pow2n_model <-  glm(BSig ~ DeltaTrend, data = neg, family = "binomial")
  
  nd <- data.frame("DeltaTrend" = seq(range(neg$DeltaTrend)[1], range(neg$DeltaTrend)[2], by = 1e-5))
  nd$P_pow1n <- predict(pow1n_model, newdata = nd, type = "response")
  nd$P_pow2n <- predict(pow2n_model, newdata = nd, type = "response")
  
  plot(nd$DeltaTrend, nd$P_pow1n, col= "blue")
  points(nd$DeltaTrend, nd$P_pow2n, col = "black")
  
}

neg <- rbind.data.frame(neg1, neg2)


length(which(tr$Sig == 0 & tr$overlaps == FALSE))
