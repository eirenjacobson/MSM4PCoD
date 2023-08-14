
# Function to initialize a population at carrying capacity
# fecundity and survival are stochastic
# returns a Z matrix of individuals alive at year 50
# last modified by EKJ 2023-08-11

initPop <- function(S0 = 0.85, S1 = 0.9, S2 = 0.95, 
                   AFR = 10, AJU = 3, ASA = 5, AMAX = 50,
                   fmax = 0.25, K = 100, nyears = 100){
  # S0 = calf survival
  # S1 = juvenile survival
  # S2 = subadult and adult survival
  # AFR = age at first reproduction
  # AJU = age at which become juvenile (subject to S1)
  # ASA = age at which become subadult (subject to S2)
  # fmax = maximum fecundity
  # K = carrying capacity
  # nyears = number of years to run to get stable age dist

  z <- 2.39 # degree of compensation
  
  # construct a vector of survival probabilities
  # note ages are zero indexed
  # i.e. calves have age 0
  S <- rep(NA, AMAX)
  S[1:AJU] <- S0
  S[(AJU+1):ASA] <- S1
  S[(ASA+1):AMAX] <- S2
  S[AMAX] <- 0
  
  # Calculate fec at equilibrium according to Eq 3 in Brandon et al
  f0 <- (1-S2) / ( (S0^AJU) * (S1^(ASA-AJU)) * (S2^(AFR-ASA)) * (1-(S2^(AMAX-AFR))) )
  
  # calculate proportion of individuals at each age class at equilibrium
  P <- rep(NA, (AMAX+1)) 
  P[1] <- 1
  
  for(i in 2:(AMAX+1)){
    P[i] <- prod(S[1:(i-1)])
  }
  
  R0 <- K/sum(P[2:length(P)])
  NPvec <- round(P * R0)
  
  # now set up a matrix of Nind x Nyears
    
  Zinit <- matrix(nrow = sum(NPvec), ncol = AMAX+1, 
                  dimnames = list(stri_rand_strings(n = sum(NPvec), length = 10)))
    
# Zinit <- matrix(nrow = sum(NPvec), ncol = AMAX+1)
  c <- 1
  for (a in 1:length(NPvec)){
    if(NPvec[a] == 0) next else{
    for (i in 1:NPvec[a]){
      Zinit[c,] <- rev(c(rep(1, a), rep(0, AMAX+1-a)))
      c <- c+1
    }}
  }
    
  return(Zinit)
    
} # end function
  