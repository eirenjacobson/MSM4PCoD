
redistributePop <- function(Na, Nb, Ka, Kb, c){
  
  # function to redistribute whales between two populations
  # Cohort based
  # Na and Nb are row vectors of population structure at some time t
  # Ka and Kb are the carrying capacities that Na and Nb are redistributing to
  # c is the degree of connectivity (0 = isolated, 1 = ideal free redistribution)
  # A is a vector of ages eligible for redistribution (typically juv/subadults)
  
  # determine how many animals would be in each cell if redistribution were perfect
  NKa <- sum(Na)/Ka
  NKb <- sum(Nb)/Kb
  
#  if(NKa == NKb){
#    # TODO: if the ratios are even, no movement happens
#  }

  
  NKt <- mean(c(NKa, NKb)) # target ratio
  
  Nma <- sum(Na)-(NKt*Ka) # how many animals need to move from a under IFR
  Nmb <- sum(Nb)-(NKt*Kb) # how many animals need to move from b under IFR
  
  if(Nma > 0){
    
    Nmac <- Nma*c # how many will actually move depends on connectivity
    Na_available <- sum(Na[A]) # how many eligible animals are available
    if(Na_available < Nmac){Na_move <- Na[A]} else { # if not enough, move them all
      Na_move <- (Na[A]/sum(Na[A])*Nmac)} # otherwise move in prop to eligible ages
    
    Na_new <- Na
    Na_new[A] <- Na[A] - Na_move
    Nb_new <- Nb
    Nb_new[A] <- Nb[A] + Na_move
    
  } else { # end if Nma > 0 
  
  if(Nmb > 0){
    
    Nmbc <- Nmb*c # how many will actually move depends on connectivity
    Nb_available <- sum(Nb[A]) # how many eligible animals are available
    if(Nb_available < Nmbc){Nb_move <- Nb[A]} else { # if not enough, move them all
      Nb_move <- (Nb[A]/sum(Nb[A])*Nmbc)} # otherwise move in prop to eligible ages
    
    Na_new <- Na
    Na_new[A] <- Na[A] + Nb_move
    Nb_new <- Nb
    Nb_new[A] <- Nb[A] - Nb_move
    
  } # end if Nma > 0
  } # end if/else
  
  Nout <- list(Na_new = Na_new, Nb_new = Nb_new)
  
  return(Nout)
}