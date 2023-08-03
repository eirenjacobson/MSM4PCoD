runIndSim <- function(nyears, cval, Ka_1, Ka_2, Kb_1, Kb_2, 
                      linetrans, linetransyrs, lt_ecv, 
                      caprecap, caprecapyrs, pcap, calfratio,
                      pam, pamyrs, pam_ecv){

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
source("./Scripts/simCalfRatio.R")

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
  
  if(calfratio == TRUE){
    if(t == caprecapyrs[1]){
      CalfDF <- data.frame()}
    if(t %in% caprecapyrs){
      calves <- simCalfRatio(Zmat = Zb_tplus1, pcap = pcap)
      adults <- nrow(filter(cr, Year == t))
      CalfDF <- rbind.data.frame(CalfDF, data.frame("Year" = t, "Calves" = calves, "Adults" = adults))
    }
  }

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

simout <- list(NSim = N, CRData = CRData, LTData = LTData, PAMData = PAMData, CalfData = CalfDF)

return(simout)

} # end function