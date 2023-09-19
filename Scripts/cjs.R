cjs <- nimbleCode({
  
  # Fixed values
  
  # Priors for parameters to be estimated by the model
  #S2 ~ dunif(0.94, 0.96)
  dS ~ dbeta(shape1 = 2, shape2 = 2)
  S2 <- 0.9 + dS*(1-0.9)
  PCap ~ dunif(0, 1)
  
  for (i in 1:Nind){
    Y[i, Find[i]:ncryears] ~ dCJS_ss(S2, PCap, len = ncryears - Find[i] + 1)
  }})