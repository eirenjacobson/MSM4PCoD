
redistributePop <- function(Za, Zb, Ka, Kb, c = 1, A = 3:10){
  
  # function to redistribute whales between two populations
  # Individual based
  # Ka and Kb are the carrying capacities that Na and Nb are redistributing to
  # c is the degree of connectivity (0 = isolated, 1 = ideal free redistribution)
  # A is a vector of ages eligible for redistribution (typically juv/subadults)

  # How many animals are in each population?
  Na <- sum(Za[,ncol(Za)])
  Nb <- sum(Zb[,ncol(Zb)])
  
  # determine how many animals would be in each cell if redistribution were perfect
  NKa <- Na/Ka
  NKb <- Nb/Kb
  
  NKt <- (Na+Nb)/(Ka+Kb) # target ratio
  
  Nma <- round((NKt*Ka) - Na) # how many animals need to move from a under IFR
  Nmb <- round((NKt*Kb) - Nb) # how many animals need to move from b under IFR

  if((Nma + Nmb) == 1 | Nma + Nmb == 0){ #if the ratios are even, no movement happens
    Za_new <- Za
    Zb_new <- Zb
  } else {
    
  # negative indicates animals need to *leave* that region
  
  if(Nma < 0 | Nmb > 0){
    
    # find out which individuals are eligible to move
    Na_alive <- which(Za[,ncol(Za)] == 1)
    # how many eligible animals are available
    Na_available <- which(rowSums(Za[Na_alive,]) %in% A)
    
    Nmac <- Nma*c # how many will actually move depends on connectivity
    Na_move <- vector() # initialize vector to hold indices of animals that move
    if(length(Na_available) <= abs(Nmac)){Na_move <- Na_available} else { # if not enough, move them all
      pmove_a <- abs(Nmac/length(Na_available)) # what is the probability of moving?
      for (i in 1:length(Na_available)){
        ifelse(rbinom(n = 1, size = 1, prob = pmove_a), Na_move <- c(Na_move, Na_available[i]), Na_move <- Na_move)
      } # end for i
    } # end else
    if(length(Na_move) == 0){
      Za_new <- Za
      Zb_new <- Zb
    } else{
    # move entries in z matrix 
    Za_new <- Za[-Na_move,]
    Zb_new <- rbind(Zb, Za[Na_move,])}
    
  } # end if Nma < 0
  
  if(Nmb < 0 | Nma > 0){
    
    # find out which individuals are eligible to move
    Nb_alive <- which(Zb[,ncol(Zb)] == 1)
    # how many eligible animals are available
    Nb_available <- which(rowSums(Zb[Nb_alive,]) %in% A)
    
    Nmbc <- Nmb*c # how many will actually move depends on connectivity
    Nb_move <- vector() # initialize vector to hold indices of animals that move
    if(length(Nb_available) <= abs(Nmbc)){Nb_move <- Nb_available} else { # if not enough, move them all
      pmove_b <- abs(Nmbc/length(Nb_available)) # what is the probability of moving?
      for (i in 1:length(Nb_available)){
        ifelse(rbinom(n = 1, size = 1, prob = pmove_b), Nb_move <- c(Nb_move, Nb_available[i]), Nb_move <- Nb_move)
      } # end for i
    } # end else
    if(length(Nb_move) == 0){
      Za_new <- Za
      Zb_new <- Zb
    } else {
    Zb_new <- Zb[-Nb_move,]
    Za_new <- rbind(Za, Zb[Nb_move,])}
    
  } # end if Nmb < 0
  
  } # end if else
  
  Nout <- list(Za_new = Za_new, Zb_new = Zb_new)
  
  return(Nout)
}