
calcPower <- function(id){
  library(dplyr)
  library(tidyr)
  library(foreach)
  library(doParallel)
  #id <- "D50_LCP_2023-06-27"
  nsim <- 100
  
  registerDoParallel(cores=10)
  #id <- "D75_LCP_MidCval_2023-07-20"
  #load(paste0("./Documents/MSM4PCoD/Results/ProcResults_", id, ".RData")) #results.out
  load(paste0("./Results/ProcResults_", id, ".RData")) #results.out
  
  ipm.trend <- data.frame()
  trend.results <- data.frame()
  power.results <- data.frame()
  
  sdf <- results.out$sdf
  ltdf <- results.out$ltdf
  pamdf <- results.out$pamdf
  rdf <- results.out$rdf
  Ndf <- results.out$Ndf
  
  tr <- foreach(i = 1:nsim, .combine = 'rbind.data.frame') %dopar% {
    for (j in unique(Ndf$Chain)){
      for (k in 1:25){ # NOTE should eventually change back to 250
        
  
        print(paste(i, j, k))
        
        psample <- dplyr::filter(Ndf, Iter == i, Chain == j, Sample == k)
        psample$Indic <- ifelse(psample$Year<=50, 0, 1)
        
        m1 <- glm(Ntot ~ Year*Indic, data = psample, family = "quasipoisson")
        ipm.trend <- rbind.data.frame(ipm.trend, 
                                          data.frame("Iter" = i, "Chain" = j, "Sample" = i, 
                                                     "Coef" = as.numeric(m1$coefficients[4])))

      } # end k (250 samples)
      } # end for j  (4 chains)
    
        a <- as.numeric(sign(quantile(ipm.trend$Coef, 0.025)) == sign(quantile(ipm.trend$Coef, 0.975)))
        trend.results <- rbind.data.frame(trend.results,
                                          data.frame("Iteration" = i, "Type" = "IPM",
                                                     "DeltaTrend" = as.numeric(quantile(ipm.trend$Coef, 0.5)),
                                                     "Sig" = a))
        ssample <- dplyr::filter(sdf, Iter == i)
        ssample$Indic <- ifelse(ssample$Year <= 50, 0, 1)*(ssample$Year - 50) 
        s_m <- glm(Ntot ~ Year + Indic, data = ssample, family = "poisson")
        trend.results <- rbind.data.frame(trend.results,
                                          data.frame("Iteration" = i, "Type" = "Sim", 
                                                     "DeltaTrend" = as.numeric(s_m$coefficients[3]),
                                                     "Sig" = as.numeric(coef(summary(s_m))[12] < 0.05)))
        
        ltsample <- dplyr::filter(ltdf, Iter == i)
        ltsample$Indic <- ifelse(ltsample$Year <= 50, 0, 1)*(ltsample$Year - 50) 
        lt_m <- glm(Nhat ~ Year + Indic, data = ltsample, family = "quasipoisson")
        trend.results <- rbind.data.frame(trend.results,
                                          data.frame("Iteration" = i, "Type" = "LT", 
                                                     "DeltaTrend" = as.numeric(lt_m$coefficients[3]),
                                                     "Sig" = as.numeric(coef(summary(lt_m))[12] < 0.05)))
        
        pamsample <- dplyr::filter(pamdf, Iter == i)
        pamsample$Indic <- ifelse(pamsample$Year <= 50, 0, 1)*(pamsample$Year-50) 
        pam_m <- glm(Nhat ~ Year + Indic, data = pamsample, family = "quasipoisson")
        trend.results <- rbind.data.frame(trend.results,
                                          data.frame("Iteration" = i, "Type" = "PAM", 
                                                     "DeltaTrend" = as.numeric(pam_m$coefficients[3]),
                                                     "Sig" = as.numeric(coef(summary(pam_m))[12] < 0.05)))
  
        trend.results
        
    } # end for.each
  
  trend.results <- tr
  
  save(trend.results, file = paste0("./Results/Trend_", id, ".RData"))
#  save(power.results, file = paste0("./Results/PowerResults", id, ".RData"))
  
} # end function

# # LT
# plot(filter(trend.results, Type == "Sim" & Period == "Post")$Trend,
#      filter(trend.results, Type == "LT" & Period == "Post")$Trend)
# 
# # PAM
# plot(filter(trend.results, Type == "Sim" & Period == "Post")$Trend,
#      filter(trend.results, Type == "PAM" & Period == "Post")$Trend)
# 
# # IPM
# plot(filter(trend.results, Type == "Sim" & Period == "Post")$Trend,
#      filter(trend.results, Type == "IPM" & Period == "Post")$Trend)
