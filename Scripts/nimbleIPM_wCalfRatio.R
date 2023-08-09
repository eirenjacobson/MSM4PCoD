
ipm <- nimbleCode({
  
  # Fixed values
  
  # Priors for parameters to be estimated by the model
  S2 ~ dunif(0.94, 0.96)
  K1_scalar ~ dunif(0, 1)
  K2_scalar ~ dunif(0, 1)
  PCap ~ dunif(0, 1)
  
  K1 <- K_lower + K1_scalar*(K_upper-K_lower)
  K2 <- K_lower + K2_scalar*(K_upper-K_lower)
  
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
  f0 <- min(fmax, (1-S2)/(S0^(AJU-1) * S1^(ASA-AJU) * S2^(AFR-ASA) * (1-S2^(AMAX-AFR-1))))
  
  # Calculate proportion in each age class under equilibrium conditions
  P[1] <- 1
  P[2] <- S0

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
  Ntot[1] <- sum(N[1, 1:AMAX])

  for (t in 1:(nyears-1)){
    ft[t] <- max(0, f0 + (fmax-f0)*(1-(sum(N[t, 2:AMAX])/K[t])^z)) # fec at t
    #N[t+1, 1] ~ dbin(ft[t], round(sum(N[t, AFR:AMAX]))) # calves
    N[t+1, 1] <- round(sum(N[t, AFR:AMAX])*ft[t])
    for (a in 2:AMAX){
      #N[t+1, a] ~ dbin(S[a-1], N[t, a-1]) # juv and adult 
      N[t+1, a] <- round(S[a-1]*N[t, a-1])
    } # end for a
    Ntot[t+1] <- sum(N[t+1, 1:AMAX])
  } # end for t
  
  # Visual survey observation model
  # including calves
  
  for (t in 1:nltyears){
     trueN[t] <- sum(N[ltyears[t],1:AMAX])
     mu_log_lt[t] <- log(trueN[t]^2/sqrt(trueN[t]^2 + ltestSD[t]^2))
     sigma_log_lt[t] <- sqrt(log(1+(ltestSD[t]^2/trueN[t]^2)))
     ltestN[t] ~ dlnorm(meanlog = mu_log_lt[t], sdlog = sigma_log_lt[t])
  }
  
  # Passive acoustic survey

  for (t in 1:npamyears){
    trueNb[t] <- sum(N[pamyears[t],1:AMAX])/2
    mu_log_pam[t] <- log(trueNb[t]^2/sqrt(trueNb[t]^2 + pamestSD[t]^2))
    sigma_log_pam[t] <- sqrt(log(1+(pamestSD[t]^2/trueNb[t]^2)))
    pamestN[t] ~ dlnorm(meanlog = mu_log_pam[t], sdlog = sigma_log_pam[t])
  }
  
  # Capture-recapture process and observation model
  
  for (i in 1:Nind){
    Y[i, Find[i]:ncryears] ~ dCJS_ss(S2, PCap, len = ncryears - Find[i] + 1)
  }
  
  for (t in 1:ncryears){
    # trials is adults observed (known), P is calf/adult ratio (at time t)
   modelCalves[t] <- sum(N[cryears[t], 1:2])
   modelAdults[t] <- sum(N[cryears[t], ASA:AMAX])
   obsAdults[t] <- sum(Y[,cryears[t]])
   NCalves[t] ~ dbinom(modelCalves[t]/modelAdults[t], obsAdults[t])
   
  }

})