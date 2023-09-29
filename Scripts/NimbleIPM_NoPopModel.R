
ipm <- nimbleCode({
  
  # Fixed values
  
  # Priors for parameters to be estimated by the model
  
  K1_scalar ~ dunif(0, 1)
  
  K1 <- K_lower + K1_scalar*(K_upper-K_lower)
  
  for (t in 1:100){
    ltestN[t] ~ dnorm(K1, sd = 50)
  }
  
  # Passive acoustic survey
  
  for (t in 1:100){
    pamestN[t] ~ dnorm(K1, sd = 50)
  }
  
  
})