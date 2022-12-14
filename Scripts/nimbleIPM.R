
ipm <- nimbleCode({
  
  # Fixed values
  
  
  # Priors for parameters to be estimated by the model
  S2 ~ dunif(0,1)
  K1 ~ dunif(175, 225)
  K2 ~ dunif(125, 175)
  # need recapture probability
  
  
  # Process model
  
  # construct a vector of survival probabilities
  #S[1:AMAX] <- rep(NA, AMAX-1)
  S[1:(AJU-1)] <- S0
  S[AJU:(ASA-1)] <- S1
  S[ASA:AMAX] <- S2
  
  # construct a vector of carrying capacity
  K[1:nyears] <- c(rep(K1, cyear - 1), rep(K2, length(nyears-cyear)))
  
  # Calculate fec at equilibrium according to Eq 3 in Brandon et al
  # Density-dependence affects fecundity according to Pella-Tomlinson
  f0 <- (1-S2)/(S0^(AJU-1)*S1^(ASA-1)*S2^(AFR-ASA-2))*(1-S2^(AMAX-AFR-2))
  
  #N[1:nyears, 1:AMAX] <- matrix(rep(NA, nyears*AMAX), ncol=AMAX) # population matrix
  
  #P <- matrix(value = NA, nrow = 1, ncol = AMAX)
  P[1] <- 1
  P[2] <- S0
  # NOTE this will cause problems if AJU is not == 3
  # could potentially add an if/else statement?
  #for (i in 3:(AJU-1)){
  #  P[i] <- S0*(S0^(i-2))
  #}
  for (i in AJU:(ASA-1)){
    P[i] <- S0^length(0:AJU)*S1^length(AJU:i)
  }
  for (i in ASA:AMAX){
    P[i] <- S0^length(0:AJU)*S1^(length(AJU:ASA)-1)*S2^((length(ASA:i)-1))
  }
  
  R0 <- K1/sum(P[2:AMAX])
  for (i in 1:AMAX){
    N[1, i] <- round(P[i]*R0)
  }

  for (t in 1:(nyears-1)){
    ft[t] <- (f0 + (fmax-f0)*(1-(sum(N[t, 2:AMAX])/K[t])^z)) # fec at t
    N[t+1, 1] ~ dbin(ft[t], round(sum(N[t, AFR:AMAX]))) # calves
    for (a in 2:AMAX){
      N[t+1, a] ~ dbin(S[a-1], N[t, a-1]) # juv and adult 
    } # end for a
  } # end for t
  
  # Visual survey observation model
  # including calves
  
  for (t in 1:nltyears){
    
    trueN[t] <- sum(N[ltyears[t],1:AMAX])
    
    ltestN[t] ~ dnorm(trueN[t], ltestSD[t])
    
  }
  
  # Capture-recapture observation model
 
  
  #for (t in cryears){
    
    
    
  #}
  
  
  
})