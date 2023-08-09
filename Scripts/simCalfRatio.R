
# Script to simulate one year of calf ratio data
# note the number of adults comes from simCapRecap
# this function just returns number of calves observed
# TODO: change this to match model formulation for obs model
simCalfRatio <- function(Zmat, pcap){
  
  # Zmat is the population matrix with nrow individuals and nyear columns
  # Pcap is the capture probability 
  ages <- rowSums(Zmat)
  alive <- which(Zmat[,ncol(Zmat)] == 1 & ages %in% 1:2)
  if(length(alive) == 0){ncalves <- 0} else {
  ncalves <- sum(rbinom(n = length(alive), size = 1, prob = pcap))}

  return(ncalves)
  
}