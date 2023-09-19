
ipm <- nimbleCode({
  
  # Fixed values
  
  # Priors for parameters to be estimated by the model

  K1_scalar ~ dunif(0, 1)

  
  S2 <- 0.95
  K1 <- K_lower + K1_scalar*(K_upper-K_lower)

  
  # Process model
  
  # construct a vector of survival probabilities
  S[1:AJU] <- S0
  S[(AJU+1):ASA] <- S1
  S[(ASA+1):(AMAX-1)] <- S2
  S[AMAX] <- 0
  
  # construct a vector of carrying capacity
  K[1:nyears] <- rep(K1, nyears)
  
  # Calculate fec at equilibrium according to Eq 3 in Brandon et al
  # Density-dependence affects fecundity according to Pella-Tomlinson
  f0 <- (1-S2) / ( (S0^AJU) * (S1^(ASA-AJU)) * (S2^(AFR-ASA)) * (1-(S2^(AMAX-AFR))) )
  # Calculate proportion in each age class at equilibrium
  
  P[1] <- 1
  
  for(i in 2:(AMAX+1)){
    P[i] <- prod(S[1:(i-1)])
  }
  
  # distribute K1 individuals according to P
  R0 <- K1/sum(P[2:(AMAX+1)])
  for (i in 1:AMAX){
    N[1, i] <- P[i]*R0
  }
  Ntot[1] <- sum(N[1, 1:AMAX])

  # run population forward in time
  for (t in 1:(nyears-1)){
    calves[t] <- N[t, 1]
    noncalves[t] <- sum(N[t, 2:(AMAX+1)])
    ft[t] <- max(0, (f0 + (fmax-f0)*(1-(sum(N[t, 2:(AMAX+1)])/K[t])^z))) # fec at t
    #N[t+1, 1] ~ dbin(ft[t], round(sum(N[t, AFR:AMAX]))) # calves
    N[t+1, 1] <- sum(N[t, (AFR+1):(AMAX+1)])*ft[t]
    for (a in 2:(AMAX+1)){
      #N[t+1, a] ~ dbin(S[a-1], N[t, a-1]) # juv and adult 
      N[t+1, a] <- S[a-1]*N[t, a-1]
    } # end for a
    Ntot[t+1] <- sum(N[t+1, 1:(AMAX+1)])
  } # end for t
  
  # Visual survey observation model
  # including calves
  
  for (t in 1:nltyears){
    trueN[t] <- Ntot[ltyears[t]]
     #mu_log_lt[t] <- log(trueN[t]^2/sqrt(trueN[t]^2 + ltestSD[t]^2))
     #sigma_log_lt[t] <- sqrt(log(1+(ltestSD[t]^2/trueN[t]^2)))
     #ltestN[t] ~ dlnorm(meanlog = mu_log_lt[t], sdlog = sigma_log_lt[t])
    #ltestN[t] ~ dlnorm(meanlog = log(trueN[t]), sdlog = log(ltestSD[t]))
    ltestN[t] ~ dnorm(trueN[t], ltestSD[t])
  }
  
  # Passive acoustic survey

  for (t in 1:npamyears){
    trueNb[t] <- Ntot[pamyears[t]]
    #mu_log_pam[t] <- log(trueNb[t]^2/sqrt(trueNb[t]^2 + pamestSD[t]))
    #sigma_log_pam[t] <- sqrt(log(1+(pamestSD[t]^2/trueNb[t]^2)))
    #pamestN[t] ~ dlnorm(meanlog = mu_log_pam[t], sdlog = sigma_log_pam[t])
    #pamestN[t] ~ dlnorm(meanlog = log(trueNb[t]), sdlog = log(pamestSD[t]))
    pamestN[t] ~ dnorm(trueNb[t], pamestSD[t])
  }


})