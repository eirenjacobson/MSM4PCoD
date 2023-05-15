
# Script to simulate one year of capture-recapture data

simCapRecap <- function(Zmat, pcap){
  
  ASA <- 5
  AMAX <- 50
  # Zmat is the population matrix with nrow individuals and nyear columns
  # Pcap is the capture probability 
  ages <- rowSums(Zmat)
  alive <- Zmat[which(Zmat[,ncol(Zmat)] == 1 & ages %in% ASA:AMAX),]
  
  IDs <- dimnames(alive[which(rbinom(n = nrow(alive), size = 1, prob = pcap) == 1),])[[1]]
  
  return(IDs)
  
}