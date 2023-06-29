
# Script to simulate one year of calf ratio data
# note the number of adults comes from simCapRecap
# this function just returns number of calves observed

simCalfRatio <- function(Zmat, pcap){
  
  # Zmat is the population matrix with nrow individuals and nyear columns
  # Pcap is the capture probability 
  ages <- rowSums(Zmat)
  alive <- Zmat[which(Zmat[,ncol(Zmat)] == 1 & ages %in% 1:2),]
  if(nrow(alive) == 0){ncalves <- 0} else {
  ncalves <- rbinom(n = nrow(alive), size = 1, prob = pcap)}

  return(ncalves)
  
}