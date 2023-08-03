
calcPower <- function(id){
  library(dplyr)
  library(tidyr)
  library(foreach)
  library(doParallel)
  #id <- "D50_LCP_2023-06-27"
  cyear <- 50
  nsim <- 100
  
  registerDoParallel(cores=10)
  
  #load(paste0("./Documents/MSM4PCoD/Results/ProcResults_", id, ".RData")) #results.out
  load(paste0("./Results/ProcResults_", id, ".RData")) #results.out
  
  trend.results <- data.frame()
  power.results <- data.frame()
  
  sdf <- results.out$sdf
  ltdf <- results.out$ltdf
  pamdf <- results.out$pamdf
  rdf <- results.out$rdf
  Ndf <- results.out$Ndf
  
  tr <- foreach(i = 1:100, .combine = 'rbind.data.frame') %dopar% {
    for (j in unique(Ndf$Chain)){
      for (k in 1:250){
        print(paste(i, j, k))
        m1 <- glm(Ntot ~ Year, 
                 data = dplyr::filter(Ndf, Year %in% 1:50, Iter == i, Chain == j, Sample == k),
                 family = "quasipoisson")
        power.results <- rbind.data.frame(power.results, 
                                          data.frame("Iter" = i, "Chain" = j, 
                                                     "Sample" = i, "Period" = "Pre",
                                                     "Coef" = as.numeric(m1$coefficients[2])))
        m2 <- glm(Ntot ~ Year, 
                 data = dplyr::filter(Ndf, Year %in% 51:100, Iter == i, Chain == j, Sample == k),
                 family = "quasipoisson")
        power.results <- rbind.data.frame(power.results, 
                                          data.frame("Iter" = i, "Chain" = j, 
                                                     "Sample" = i, "Period" = "Post",
                                                     "Coef" = as.numeric(m2$coefficients[2])))
      } #end k
      } # end for j and k (4 chains x 250 samples each)
    
    save(power.results, file = paste0("./Results/PowerResults_", id, ".RData"))
    
  
        precoefs <- dplyr::filter(power.results, Period == "Pre", Iter == i)$Coef
        postcoefs <- dplyr::filter(power.results, Period == "Post", Iter == i)$Coef
        a1 <- as.numeric(sign(quantile(precoefs, 0.025)) == sign(quantile(precoefs, 0.975)))
        a2 <- as.numeric(sign(quantile(postcoefs, 0.025)) == sign(quantile(postcoefs, 0.975)))
        trend.results <- rbind.data.frame(trend.results,
                                          data.frame("Iteration" = i, "Type" = "Model",
                                                     "Period" = "Pre", 
                                                     "Trend" = as.numeric(quantile(precoefs, 0.5)),
                                                     "Sig" = a1),
                                          data.frame("Iteration" = i, "Type" = "Model",
                                                     "Period" = "Post",
                                                     "Trend" = as.numeric(quantile(postcoefs, 0.5)),
                                                     "Sig" = a2))
  
    
    s_pre <- dplyr::filter(sdf, Iter == i, Year <= cyear)
    s_pre_m <- glm(Ntot ~ Year, data = s_pre, family = "poisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "Sim", "Period" = "Pre",
                                                 "Trend" = as.numeric(s_pre_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(s_pre_m))[8] < 0.05)))
  
    s_post <- dplyr::filter(sdf, Iter == i, Year > cyear)
    s_post_m <- glm(Ntot ~ Year, data = s_post, family = "poisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "Sim", "Period" = "Post",
                                                 "Trend" = as.numeric(s_post_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(s_post_m))[8] < 0.05)))
  
    lt_pre <- dplyr::filter(ltdf, Iter == i, Year <= cyear)
    lt_pre_m <- glm(Nhat ~ Year, data = lt_pre, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "LT", "Period" = "Pre",
                                                 "Trend" = as.numeric(lt_pre_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(lt_pre_m))[8] < 0.05)))
  
    lt_post <- dplyr::filter(ltdf, Iter == i, Year > cyear)
    lt_post_m <- glm(Nhat ~ Year, data = lt_post, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "LT", "Period" = "Post",
                                                 "Trend" = as.numeric(lt_post_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(lt_post_m))[8] < 0.05)))
  
    pam_pre <- dplyr::filter(ltdf, Iter == i, Year <= cyear)
    pam_pre_m <- glm(Nhat ~ Year, data = pam_pre, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "PAM", "Period" = "Pre",
                                                 "Trend" = as.numeric(pam_pre_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(pam_pre_m))[8] < 0.05)))
  
    pam_post <- dplyr::filter(pamdf, Iter == i, Year > cyear)
    pam_post_m <- glm(Nhat ~ Year, data = pam_post, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "PAM", "Period" = "Post",
                                                 "Trend" = as.numeric(pam_post_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(pam_post_m))[8] < 0.05)))
  
    ipm_pre <- dplyr::filter(rdf, Iter == i, Year <= cyear)
    ipm_pre_m <- glm(Ntot ~ Year, data = ipm_pre, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "IPM", "Period" = "Pre",
                                                 "Trend" = as.numeric(ipm_pre_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(ipm_pre_m))[8] < 0.05)))
    
    ipm_post <- dplyr::filter(rdf, Iter == i, Year > cyear)
    ipm_post_m <- glm(Ntot ~ Year, data = ipm_post, family = "quasipoisson")
    trend.results <- rbind.data.frame(trend.results,
                                      data.frame("Iteration" = i,
                                                 "Type" = "IPM", "Period" = "Post",
                                                 "Trend" = as.numeric(ipm_post_m$coefficients[2]),
                                                 "Sig" = as.numeric(coef(summary(ipm_post_m))[8] < 0.05)))
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
