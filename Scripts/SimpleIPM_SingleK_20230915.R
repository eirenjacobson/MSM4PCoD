library(nimble)
library(runjags)
library(ggplot2)


source("./Scripts/simPop.R")
source("./Scripts/simSurvey.R")
source("./Scripts/nimbleIPM_Simple_SingleK_20230918.R")

targetN <- 500
nyears <- 100
ltcv <- 0.1
pamcv <- 0.1

set.seed(20230920)

results <- list()

for (i in 1:10){

pop <- simPop(K = targetN, nyears = nyears)

N.true <- rowSums(pop$N)

LT.data <- vector()
PAM.data <- vector()
for (j in 1:100){
  LT.data[j] <- simSurvey(N.true[j], ltcv)
  PAM.data[j] <- simSurvey(N.true[j], pamcv)
}

Ndefault <- simPop(K = targetN, nyears = nyears)  

nimbleData <- list(ltestN = LT.data, ltestSD = LT.data*ltcv, 
                   pamestN = PAM.data, pamestSD = PAM.data*pamcv)

nimbleConstants <- list(nyears = nyears, S0 = 0.85, S1 = 0.9, AFR = 10,
                        AJU = 3, ASA = 5, AMAX = 50, fmax = 0.2,  
                        z = 2.39, 
                        ltyears = 1:nyears, nltyears = length(1:nyears),
                        pamyears = 1:nyears, npamyears = length(1:nyears),
                        K_lower = 400,
                        K_upper = 600)

nimbleInits <- list(K1_scalar = 1, 
                    N = round(Ndefault$N), Ntot = round(rowSums(round(Ndefault$N))))

nimbleParams <- list("K1",  "noncalves", "calves")


model <- nimbleModel(code = ipm,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        constants = nimbleConstants, data = nimbleData,
                        thin = 10, niter = 50000, nburnin = 40000, nchains = 4,
                        samplesAsCodaMCMC = TRUE)

results[[i]] <- nimbleOut

}

K.results <- data.frame("i" = 1:10, "K1_Med" = NA, "K1_LCI" = NA, "K1_UCI" = NA)
N.results <- data.frame()

for (i in 1:10){
  
  single.mcmc <- combine.mcmc(results[[i]])
  
  K.results[i,] <- c(i, summary(single.mcmc)$quantiles[1,3], 
                     summary(single.mcmc)$quantiles[1,1],
                     summary(single.mcmc)$quantiles[1,5])
  N.results <- rbind.data.frame(N.results, data.frame("i" = i, "Year" = 1:(nyears-1), 
                          "calves" = summary(single.mcmc)$quantiles[2:100,3],
                          "noncalves" = summary(single.mcmc)$quantiles[101:199, 3]))
  
}

ggplot(data = N.results)+
  geom_line(aes(x=Year, y = noncalves, group = i)) +
  geom_line(aes(x=Year, y = noncalves+calves, group = i))+
  #geom_line(data = K.results, aes(x=1:1, y = K1_Med))+
#  geom_segment(data = K.results, aes(x = 0, xend = 100, y = K.results$K1, yend = K.results$K1), color = "blue")+
  facet_wrap(~i)+
  theme_bw()



K1_svals <- seq(0.4, 0.6, length.out = 100)

ll <- rep(NA, 100)

for (i in 1:100){
  
  model$setInits(list("K1_scalar" = K1_svals[i]))
  ll[i] <- calculate(model)
  
}

plot(K1_svals, ll)
abline(v = K1_svals[which.max(ll)])
abline(v=0.5, col = "red")
