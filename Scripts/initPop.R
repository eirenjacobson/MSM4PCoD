
# Function to initialize a population at carrying capacity
# fecundity and survival are stochastic
# returns a Z matrix of individuals alive at year 50
# last modified by EKJ 2022-11-17

initPop <- function(S0 = 0.8, S1 = 0.85, S2 = 0.95, 
                   AFR = 10, AJU = 3, ASA = 5, AMAX = 50,
                   fmax = 0.2, K = 100, nyears = 100){
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
  S <- rep(NA, AMAX-1)
  S[1:(AJU-1)] <- S0
  S[AJU:(ASA-1)] <- S1
  S[ASA:AMAX] <- S2
  
  # Calculate fec at equilibrium according to Eq 3 in Brandon et al
  # Density-dependence affects fecundity according to Pella-Tomlinson
  f0 <- (1-S2)/(S0^(AJU-1)*S1^(ASA-1)*S2^(AFR-ASA-2))*(1-S2^(AMAX-AFR-2))
    
    N <- matrix(rep(NA, nyears*AMAX), ncol=AMAX) # population matrix
    
    P <- rep(NA, AMAX) 
    P[1] <- 1
    P[2] <- S0
    for (i in 3:(AJU-1)){
      P[i] <- S0*(S0^(i-2))
    }
    for (i in AJU:(ASA-1)){
      P[i] <- S0^length(0:AJU)*S1^length(AJU:i)
    }
    for (i in ASA:AMAX){
      P[i] <- S0^length(0:AJU)*S1^(length(AJU:ASA)-1)*S2^((length(ASA:i)-1))
    }
    
    R0 <- K/sum(P[2:length(P)])
    N[1, ] <- round(P * R0)
    
    for (t in 1:(nyears-1)){
      ft <- (f0 + (fmax-f0)*(1-(sum(N[t, 2:AMAX])/K)^z)) # fec at t
      N[t+1, 1] <- rbinom(n = 1, size = round(sum(N[t, AFR:AMAX])), prob = ft) # calves
      for (a in 2:AMAX){
        N[t+1, a] <- rbinom(n = 1, size = N[t, a-1], prob = S[a-1]) # juv and adult 
      } # end for a
    } # end for t
    
    # now set up a matrix of Nind x Nyears
    
    Zinit <- matrix(nrow = sum(N[nyears,]), ncol = max(which(N[nyears,] > 0)),
                   dimnames = list(stri_rand_strings(n = sum(N[nyears,]), length = 10)))
    
    c <- 1
    for (a in 1:ncol(Zinit)){
      if(N[nyears,a] == 0){next} else{
      for (i in 1:N[nyears,a]){
        Zinit[c, 1:(ncol(Zinit)-a)] <- 0
        Zinit[c, (ncol(Zinit)-a+1):ncol(Zinit)] <- 1
        c <- c+1
      } 
      }}
    
    # in case no ind alive at year 50 are 50 years old
    # add some padding zeroes to the beginning of the matrix
    # to maintain compatibility with other pops
    if(ncol(Zinit < 51)){
      ac <- 51 - ncol(Zinit) 
      Zinit <- cbind(matrix(rep(0, nrow(Zinit)*ac), ncol = ac), Zinit)
      
    }
    
    return(Zinit)
    
} # end function
  