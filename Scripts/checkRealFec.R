

# Set vector of carrying capacities
Ka <- c(rep(Kb_1, nyears/2), rep(Kb_2, nyears/2))
Kb <- c(rep(Kb_1, nyears/2), rep(Kb_2, nyears/2))

# Set counter for number of animals in each area
Na_t <- rep(NA, nyears)
Nb_t <- rep(NA, nyears)

nadults <- rep(NA, nyears)
ncalves <- rep(NA, nyears)

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
  
  nadults[t] <- length(which(rowSums(Za_tplus1) > AFR & Za_tplus1[,ncol(Za_tplus1)] == 1))
  ncalves[t] <- length(which(rowSums(Za_tplus1) == 1 & Za_tplus1[,ncol(Za_tplus1)] == 1))
  
} # end for t
