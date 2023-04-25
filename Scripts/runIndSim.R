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
source("./Scripts/simLineTransect.R")

############################

nyears <- 100
set.c <- 1
caprecapyrs <- 40:60
linetransyrs <- seq(2,99, by = 1)

Area_A <- 2500
Area_B <- 2500

LT_A <- 50
LT_B <- 50

ESW <- 2

############################

set.seed(20221116)

# Set vector of carrying capacities
Ka <- rep(100, nyears)
#Kb <- c(rep(100, nyears/2), rep(50, nyears/2))
Kb <- Ka
#plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(-10, 10))

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
  
  Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t], c = set.c)
  # These z-matrices will be the starting point for the next iteration
  Za_tplus1 <- Z_new$Za_new
  Zb_tplus1 <- Z_new$Zb_new
  
  # Fill in vectors with total number of animals in each area
  Na_t[t] <- sum(Za_tplus1[,ncol(Za_tplus1)])
  Nb_t[t] <- sum(Zb_tplus1[,ncol(Zb_tplus1)])
  
  # now simulate surveys
  if (t %in% caprecapyrs){

    if(t == caprecapyrs[1]){
      caprecap <- data.frame("ID" =simCapRecap(Zmat = Zb_tplus1, Pcap = 0.1),
                             "Year" = t, "Cap" = 1)
    } else {
    caprecap <- rbind.data.frame(caprecap, 
                                 data.frame("ID" =simCapRecap(Zmat = Zb_tplus1, Pcap = 0.1),
                                            "Year" = t, "Cap" = 1))}
  } # end caprecap

  if (t %in% linetransyrs){
    if(t == linetransyrs[1]){
      
      Nhat <- data.frame("Year" = linetransyrs, "Nhat" = NA, "SD" = NA)
      LTEst <- simLineTransect(N = Na_t[t] + Nb_t[t], L = LT_A + LT_B, A = Area_A + Area_B)
      Nhat$Nhat[linetransyrs == t] <- LTEst$Nhat
      Nhat$SD[linetransyrs == t] <- LTEst$SD} else {

    LTEst <- simLineTransect(N = Na_t[t] + Nb_t[t], L = LT_A + LT_B, A = Area_A + Area_B)
    Nhat$Nhat[linetransyrs == t] <- LTEst$Nhat
    Nhat$SD[linetransyrs == t] <- LTEst$SD
  }} # end line transect
  
} # end for t

caphist <- caprecap %>% 
  pivot_wider(id_cols = ID, names_from = Year, values_from = Cap, values_fill = 0)

# Process caphist for nimbleIPM

F <- caprecap %>% 
  group_by(ID) %>% 
  filter(Year == min(Year)) %>% 
  ungroup() %>%
  select(Year) %>% as.vector()

Y <- as.matrix(caphist[,2:ncol(caphist)])
# take out individuals only seen on last occasion
Y <- Y[-which(rowSums(Y) == 1 & Y[,21] == 1),]

Zmat <- matrix(rep(0, length(F)*ncol(Y)))

# Note these are the tplus1 matrices, not the Na/Nb matrices
# N <- data.frame("Year" = c(1:ncol(Za_tplus1), 1:ncol(Zb_tplus1)), 
#                 "N" = c(colSums(Za_tplus1), colSums(Zb_tplus1)), 
#                 "Region" = c(rep("A", ncol(Za_tplus1)), rep("B", ncol(Zb_tplus1))),
#                 "K" = c(rep(100, 51), Ka, rep(100, 51), Kb))

N <- data.frame("Year" = rep(1:nyears, 2),
                "N" = c(Na_t, Nb_t), 
                "Region" = c(rep("A", length(Na_t)), rep("B", length(Nb_t))),
                "K" = c(Ka, Kb))

ggplot(N) +
  geom_line(aes(x=Year, y = N, color = Region)) +
  geom_line(aes(x=Year, y = K, color = Region), linetype = "dashed") +
 # xlim(50, 150) +
  theme_bw() 

N %>% group_by(Year) %>% summarize(Ntot = sum(N)) %>%
  ggplot()+
  geom_line(aes(x=Year, y = Ntot))+
  geom_point(data = Nhat, aes(x=Year, y = Nhat))


###################

source("./Scripts/nimbleIPM.R")
source("./Scripts/simPop.R")

Ndefault <- simPop(K = 225, new = TRUE, nyears = 100)

nimbleData <- list(ltestN = Nhat$Nhat, ltestSD = rep(10, length(Nhat$Year)), 
                   Y = Y)

nimbleConstants <- list(cyear = 50, nyears = 100, S0 = 0.8, S1 = 0.85, AFR = 10,
                        AJU = 3, ASA = 5, AMAX = 50, fmax = 0.2, nyears = 100, 
                        z = 2.39, nltyears = length(Nhat$Year), ltyears = Nhat$Year,
                        PCap = 0.1, ncryears = length(caprecapyrs), Nind = nrow(Y),
                        Find = F$Year - F$Year[1] + 1)

nimbleInits <- list(S2 = 0.95, K1 = 225, K2 = 225, N = round(Ndefault))

nimbleParams <- list("S2", "K1", "K2", "trueN")

set.seed(20220401) 

model <- nimbleModel(code = ipm,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, #constants = nimbleCi, inits = nimbleInits,
                        monitors = nimbleParams, 
                        thin = 10, niter = 1000, nburnin = 500, nchains = 4)

