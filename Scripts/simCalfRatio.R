
# Script to simulate one year of calf ratio data
# note the number of adults comes from simCapRecap
# this function just returns number of calves observed

simCalfRatio <- function(Zmat, nadults){
  
  ASA <- 5
  AMAX <- 50
  
  # Zmat is the population matrix with nrow individuals and nyear columns
  # Pcap is the capture probability 
  ages <- rowSums(Zmat)
  calves <- length(which(Zmat[,ncol(Zmat)] == 1 & ages %in% 1:2))
  noncalves <- length(which(Zmat[,ncol(Zmat)] == 1 & ages %in% ASA:AMAX))
  if(length(calves) == 0){obs <- 0} else {
  obs <- rbinom(prob = calves/noncalves, n = 1, size = nadults)
  }
  return(obs)
  
}