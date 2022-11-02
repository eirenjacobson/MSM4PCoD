
# Script to simulate one year of capture-recapture data

simCapRecap <- function(Zmat, Pcap){
  
  # Zmat is the population matrix with nrow individuals and nyear columns
  # Pcap is the capture probability

  alive <- Zmat[which(Zmat[,ncol(Zmat)] == 1),]
  
  IDs <- dimnames(alive[which(rbinom(n = nrow(alive), size = 1, prob = Pcap) == 1),])[[1]]
  
  return(IDs)
  
}