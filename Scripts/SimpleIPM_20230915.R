library(nimble)
library(runjags)
library(ggplot2)
library(dplyr)
library(tidyr)

source("./Scripts/simPop_Kvec.R")
source("./Scripts/simSurvey.R")
source("./Scripts/nimbleIPM_Simple_20230915.R")

#targetN <- 500
nyears <- 100
ltcv <- 0.1
pamcv <- 0.1

set.seed(20230915)

power.results <- list()

Kvectors <- list()
Kvectors[[1]] <- c(rep(500, 50), rep(550, 50))
Kvectors[[2]] <- c(rep(500, 50), rep(540, 50))
Kvectors[[3]] <- c(rep(500, 50), rep(530, 50))
Kvectors[[4]] <- c(rep(500, 50), rep(520, 50))
Kvectors[[5]] <- c(rep(500, 50), rep(510, 50))
Kvectors[[6]] <- rep(500, 100)
Kvectors[[7]] <- c(rep(500, 50), rep(490, 50))
Kvectors[[8]] <- c(rep(500, 50), rep(480, 50))
Kvectors[[9]] <- c(rep(500, 50), rep(470, 50))
Kvectors[[10]] <- c(rep(500, 50), rep(460, 50))
Kvectors[[11]] <- c(rep(500, 50), rep(450, 50))




for (k in 1:11){
for (i in 1:10){

pop <- simPop(nyears = nyears, K = Kvectors[[k]])

N.true <- rowSums(pop$N)

LT.data <- vector()
PAM.data <- vector()
for (j in 1:100){
  LT.data[j] <- simSurvey(N.true[j], ltcv)
  PAM.data[j] <- simSurvey(N.true[j], pamcv)
}

Ndefault <- simPop(K = Kvectors[[k]], nyears = nyears)  

nimbleData <- list(ltestN = LT.data, ltestSD = LT.data*ltcv, 
                   pamestN = PAM.data, pamestSD = PAM.data*pamcv)

nimbleConstants <- list(cyear = 50, nyears = nyears, S0 = 0.85, S1 = 0.9, AFR = 10,
                        AJU = 3, ASA = 5, AMAX = 50, fmax = 0.2,  
                        z = 2.39, 
                        ltyears = 1:nyears, nltyears = length(1:nyears),
                        pamyears = 1:nyears, npamyears = length(1:nyears),
                        K_lower = 50,
                        K_upper = 1050)

nimbleInits <- list(K1_scalar = 0.5, 
                    K2_scalar = 0.5, 
                    N = round(Ndefault$N), Ntot = round(rowSums(round(Ndefault$N))))

nimbleParams <- list("K1", "K2", "Ntot")



model <- nimbleModel(code = ipm,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        constants = nimbleConstants, data = nimbleData,
                        thin = 10, niter = 50000, nburnin = 40000, nchains = 1,
                        samplesAsCodaMCMC = TRUE)

results[[i]] <- nimbleOut

}

K.results <- data.frame("iter" = 1:10, "K1" = NA, "K1_LCI" = NA, "K1_UCI" = NA, 
                        "K2" = NA, "K2_LCI" = NA, "K2_UCI" = NA )
N.results <- data.frame()

for (i in 1:10){
  
  single.mcmc <- combine.mcmc(results[[i]])
  
  K.results[i,] <- c(i, summary(single.mcmc)$quantiles[1,3], 
                     summary(single.mcmc)$quantiles[1,1],
                     summary(single.mcmc)$quantiles[1,5],
                    summary(single.mcmc)$quantiles[2,3],
                    summary(single.mcmc)$quantiles[2,1],
                    summary(single.mcmc)$quantiles[2,5])
  N.results <- rbind.data.frame(N.results, data.frame("iter" = i, "Year" = 1:(nyears-1), 
                          "Ntot" = summary(single.mcmc)$quantiles[3:101,3]))
  
}

# ggplot(K.results) +
#   geom_errorbar(aes(x = iter, ymin = K1_LCI, ymax = K1_UCI), color = "blue") +
#   geom_errorbar(aes(x=iter, ymin = K2_LCI, ymax = K2_UCI), color = "red") +
#   geom_hline(aes(yintercept = 500))+
#  # ylim(c(0, NULL))+
#   theme_bw()+
#   xlim(c(0,11))
# 
# ggplot()+
#   geom_line(data = N.results, aes(x=Year, y = Ntot, group = iter)) +
# #  geom_segment(aes(x = 0, xend = 50, y = K.results$K1, yend = K.results$K1), color = "blue") +
# #  geom_segment(aes(x=50, xend=100, y = K.results$K2, yend = K.results$K2), color = "red")+
# #  facet_grid(~iter) +
#   theme_bw()


# Checking power 

trend.results <- data.frame()
ipm.trend <- data.frame()

for (i in 1:10){
  
  Ndf <- as.data.frame(results[[i]]) %>%
    mutate(Sample = 1:nrow(results[[i]])) %>%
    select(Sample, which(substr(names(as.data.frame(results[[i]])), 1, 4) == "Ntot")) %>%
    rename_with(~gsub("Ntot", "", .x, fixed = TRUE))  %>%
    rename_with(~gsub("[", "", .x, fixed = TRUE)) %>%
    rename_with(~gsub("]", "", .x, fixed = TRUE)) %>%
    pivot_longer(cols = 2:101, names_to = "Year", values_to = "Ntot") %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate("Iter" = i)
  
  for (j in 1:1000){
    
    psample <- dplyr::filter(Ndf, Iter == i, Sample == j)
    psample$Indic <- ifelse(psample$Year<=50, 0, 1)
    
    m1 <- glm(Ntot ~ Year*Indic, data = psample, family = "quasipoisson")
    ipm.trend <- rbind.data.frame(ipm.trend, 
                                  data.frame("Iter" = i, "Sample" = j, 
                                             "Coef" = as.numeric(m1$coefficients[4])))
    
} # end for j  (1000 samples)
  
  a <- as.numeric(sign(quantile(ipm.trend$Coef, 0.025)) == sign(quantile(ipm.trend$Coef, 0.975)))
  trend.results <- rbind.data.frame(trend.results,
                                    data.frame("Iteration" = i, "Type" = "IPM",
                                               "DeltaTrend" = as.numeric(quantile(ipm.trend$Coef, 0.5)),
                                               "Sig" = a))
  
} # end for i

power.results[[k]] <- list(trend.results = trend.results, K.results = K.results, N.results = N.results)

} # end for k

save(power.results, file = "./Data/powerresults_SimpleIPM_20230929.RData")

megatrend <- data.frame()
for (i in 1:5){
  
  r <- power.results[[i]]$trend.results
  r$simchange <- (Kvectors[[i]][100]-500)/500
  megatrend <- rbind.data.frame(megatrend, r)
  
}

megatrend %>% group_by(simchange) %>% summarize(power = sum(Sig)/10)

# for (i in 1:11){
# ggplot(power.results[[11]]$K.results) +
#   geom_errorbar(aes(x = iter, ymin = K1_LCI, ymax = K1_UCI), color = "blue") +
#   geom_errorbar(aes(x=iter, ymin = K2_LCI, ymax = K2_UCI), color = "red") +
#   geom_hline(aes(yintercept = 500))+
#  # ylim(c(0, NULL))+
#   theme_bw()+
#   xlim(c(0,11))
# }