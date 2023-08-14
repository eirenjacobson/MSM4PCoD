runIndSim <- function(nyears, cval, Ka_1, Ka_2, Kb_1, Kb_2, 
                      linetrans, linetransyrs, lt_ecv, 
                      caprecap, caprecapyrs, pcap,
                      pam, pamyrs, pam_ecv, pars){

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
Ka <- c(rep(Ka_1, pars$cyear), rep(Ka_2, pars$nyears-pars$cyear))
Kb <- c(rep(Kb_1, pars$cyear), rep(Kb_2, pars$nyears-pars$cyear))

# Set counter for number of animals in each area
Na_t <- rep(NA, nyears)
Nb_t <- rep(NA, nyears)

for (t in 1:nyears){
  # if it's the first year, initialize the population
  if (t == 1){
    Za_t <- projPop(Zinit = initPop(), nyears = 1, K = Ka[t])
    Zb_t <- projPop(Zinit = initPop(), nyears = 1, K = Kb[t])
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
  if(caprecap == TRUE){
  if (t %in% caprecapyrs){
    if(t == caprecapyrs[1]){
      IDs <- simCapRecap(Zmat = Zb_tplus1, pcap = pcap)
      if (length(IDs) == 0){cr <- data.frame()} else {
      cr <- data.frame("ID" = IDs, "Year" = t, "Cap" = 1)}
    } else {
      IDs <- simCapRecap(Zmat = Zb_tplus1, pcap = pcap)
      if (length(IDs) == 0){ cr <- cr} else{
      cr <- rbind.data.frame(cr, 
                                  data.frame("ID" = IDs, "Year" = t, "Cap" = 1))}}
  }} # end caprecap

  if(linetrans == TRUE){
  if(t == linetransyrs[1]){
    {LTNhat <- data.frame("Year" = linetransyrs, "Nhat" = NA, "SD" = NA)}}
  if (t %in% linetransyrs){
        LTEst <- simSurvey(N = Na_t[t] + Nb_t[t], CV = lt_ecv)
        LTNhat$Nhat[linetransyrs == t] <- LTEst
        LTNhat$SD[linetransyrs == t] <- lt_ecv*LTEst
  }} # end line transect
  
  if(pam == TRUE){
  if(t == pamyrs[1]){
    {PAMNhat <- data.frame("Year" = pamyrs, "Nhat" = NA, "SD" = NA)}}
  if (t %in% pamyrs){
    PAMEst <- simSurvey(N = Nb_t[t], CV = pam_ecv)
    PAMNhat$Nhat[pamyrs == t] <- PAMEst
    PAMNhat$SD[pamyrs == t] <- pam_ecv*PAMEst
  }} # end passive acoustic monitoring
  
} # end for t

N <- data.frame("Year" = rep(1:nyears, 2),
                "N" = c(Na_t, Nb_t), 
                "Region" = c(rep("A", length(Na_t)), rep("B", length(Nb_t))),
                "K" = c(Ka, Kb))

if(caprecap == TRUE){CRData <- cr}else{CRData <- NULL}
if(linetrans == TRUE){LTData <- LTNhat}else{LTData <- NULL}
if(pam == TRUE){PAMData <- PAMNhat}else{PAMData <- NULL}

simout <- list(NSim = N, CRData = CRData, LTData = LTData, PAMData = PAMData)

return(simout)

} # end function