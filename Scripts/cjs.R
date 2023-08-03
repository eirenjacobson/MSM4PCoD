cjs <- nimbleCode({
  
  # Fixed values
  
  # Priors for parameters to be estimated by the model
  S2 ~ dunif(0.94, 0.96)
  PCap ~ dunif(0, 1)
  
  for (i in 1:Nind){
    Y[i, Find[i]:ncryears] ~ dCJS_ss(S2, PCap, len = ncryears - Find[i] + 1)
  }})