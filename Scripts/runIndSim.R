runIndSim <- function(nyears, cval, Ka_1, Ka_2, Kb_1, Kb_2, 
                      linetransyrs, lt_ecv, 
                      caprecapyrs, pcap,
                      nchains, thin, niter, nburnin){

####
library(dplyr)
library(tidyr)
library(ggplot2)
library(nimble)
library(nimbleEcology)
####

source("./Scripts/initPop.R")
source("./Scripts/projPop.R")
source("./Scripts/redistributePopInd.R")
source("./Scripts/simCapRecap.R")
source("./Scripts/simSurvey.R")

############################
# 
# # years to run simulation
# nyears <- 100
# # connectivity parameter
# cval <- 1
# # years in which capture-recapture surveys happen
# caprecapyrs <- 40:60
# # capture probability 
# pcap <- 0.2
# # years in which line-transect surveys happen
# linetransyrs <- 40:60
# ecv <- 0.2


############################

# Set vector of carrying capacities
Ka <- c(rep(Kb_1, nyears/2), rep(Kb_2, nyears/2))
Kb <- c(rep(Kb_1, nyears/2), rep(Kb_2, nyears/2))

# Set counter for number of animals in each area
Na_t <- rep(NA, nyears)
Nb_t <- rep(NA, nyears)

for (t in 1:nyears){
  # if it's the first year, initialize the population
  if (t == 1){
    Za_t <- projPop(Zinit = initPop(), nyears = 1)
    Zb_t <- projPop(Zinit = initPop(), nyears = 1)
  } else {
    Za_t <- projPop(Zinit = Za_tplus1, nyears = 1, K = Ka[t])
    Zb_t <- projPop(Zinit = Zb_tplus1, nyears = 1, K = Kb[t])
  }
  
  Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t], c = cval)
  # These z-matrices will be the starting point for the next iteration
  Za_tplus1 <- Z_new$Za_new
  Zb_tplus1 <- Z_new$Zb_new
  
  # Fill in vectors with total number of animals in each area
  Na_t[t] <- sum(Za_tplus1[,ncol(Za_tplus1)])
  Nb_t[t] <- sum(Zb_tplus1[,ncol(Zb_tplus1)])
  
  # now simulate surveys
  if (t %in% caprecapyrs){

    if(t == caprecapyrs[1]){
      IDs <- simCapRecap(Zmat = Zb_tplus1, pcap = pcap)
      if (length(IDs) == 0){caprecap <- data.frame()} else {
      caprecap <- data.frame("ID" = IDs, "Year" = t, "Cap" = 1)}
    } else {
      IDs <- simCapRecap(Zmat = Zb_tplus1, pcap = pcap)
      if (length(IDs) == 0){ caprecap <- caprecap} else{
      caprecap <- rbind.data.frame(caprecap, 
                                  data.frame("ID" = IDs, "Year" = t, "Cap" = 1))}}
  } # end caprecap

  if(t == linetransyrs[1]){
    
    {Nhat <- data.frame("Year" = linetransyrs, "Nhat" = NA, "SD" = NA)}}
  if (t %in% linetransyrs){
        LTEst <- simSurvey(N = Na_t[t] + Nb_t[t], CV = lt_ecv)
        Nhat$Nhat[linetransyrs == t] <- LTEst
        Nhat$SD[linetransyrs == t] <- lt_ecv*LTEst
  } # end line transect
  
} # end for t

caphist <- complete(caprecap, ID, Year = caprecapyrs, fill = list(Cap = 0)) %>% 
  pivot_wider(id_cols = ID, names_from = Year, values_from = Cap)

# Process caphist for nimbleIPM


Y <- as.matrix(caphist[,2:ncol(caphist)])
# take out individuals only seen on last occasion
Y <- Y[-which(rowSums(Y) == 1 & Y[,ncol(Y)] == 1),]

F <- rep(0, nrow(Y))
for (i in 1:nrow(Y)){
  F[i] <- min(which(Y[i,] == 1))
}

#Zmat <- matrix(rep(0, length(F)*ncol(Y)))

# Note these are the tplus1 matrices, not the Na/Nb matrices
# N <- data.frame("Year" = c(1:ncol(Za_tplus1), 1:ncol(Zb_tplus1)), 
#                 "N" = c(colSums(Za_tplus1), colSums(Zb_tplus1)), 
#                 "Region" = c(rep("A", ncol(Za_tplus1)), rep("B", ncol(Zb_tplus1))),
#                 "K" = c(rep(100, 51), Ka, rep(100, 51), Kb))

N <- data.frame("Year" = rep(1:nyears, 2),
                "N" = c(Na_t, Nb_t), 
                "Region" = c(rep("A", length(Na_t)), rep("B", length(Nb_t))),
                "K" = c(Ka, Kb))


source("./Scripts/nimbleIPM.R")
source("./Scripts/simPop.R")

Ndefault <- simPop(K = round(max(Nhat$Nhat)), new = TRUE, nyears = nyears)

nimbleData <- list(ltestN = Nhat$Nhat, ltestSD = rep(10, length(Nhat$Year)), 
                   Y = Y)

nimbleConstants <- list(cyear = 50, nyears = nyears, S0 = 0.8, S1 = 0.85, AFR = 10,
                        AJU = 3, ASA = 5, AMAX = 50, fmax = 0.25, nyears = 100, 
                        z = 2.39, nltyears = length(Nhat$Year), ltyears = Nhat$Year,
                        pcap = pcap, ncryears = length(caprecapyrs), Nind = nrow(Y),
                        Find = F,
                        K1_upper = round(max(Nhat$Nhat)*2), 
                        K2_upper = round(max(Nhat$Nhat)*2))

nimbleInits <- list(S2 = 0.95, 
                    K1 = round(max(Nhat$Nhat)), 
                    K2 = round(max(Nhat$Nhat)), 
                    N = round(Ndefault), Ntot = rowSums(Ndefault))

nimbleParams <- list("S2", "K1", "K2", "Ntot", "ft", "f0")

model <- nimbleModel(code = ipm,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        thin = thin, niter = niter, nburnin = nburnin, nchains = nchains,
                        samplesAsCodaMCMC = TRUE)

out <- list(N = N, nimbleOut = nimbleOut, ltestN = Nhat$Nhat)

return(out)


} # end function