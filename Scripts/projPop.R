
# Function to project Ziphiid-ish population and track individuals over time
# Operates by INDIVIDUALS
# Last modified by EKJ 2022-08-23

projPop <- function(S0 = 0.8, S1 = 0.85, S2 = 0.95, 
                   AFR = 10, AJU = 3, ASA = 5, AMAX = 50,
                   fmax = 0.2, K = 1000, nyears = 1, Zinit = NA){
  
  
  library(stringi)
  # S0 = calf survival
  # S1 = juvenile survival
  # S2 = subadult and adult survival
  # AFR = age at first reproduction
  # AJU = age at which become juvenile (subject to S1)
  # ASA = age at which become subadult (subject to S2)
  # fmax = maximum fecundity
  # K = carrying capacity
  # Zinit is the output of initPop or of this function + redistribution
  
  z <- 2.39 # degree of compensation
  
  # construct a vector of survival probabilities
  S <- rep(NA, AMAX-1)
  S[1:(AJU-1)] <- S0
  S[AJU:(ASA-1)] <- S1
  S[ASA:AMAX] <- S2
  
  # Calculate fec at equilibrium according to Eq 3 in Brandon et al
  # Density-dependence affects fecundity according to Pella-Tomlinson
  f0 <- (1-S2)/(S0^(AJU-1)*S1^(ASA-1)*S2^(AFR-ASA-2))*(1-S2^(AMAX-AFR-2))
  
  # generate a Z matrix with the initial individuals in the upper left
  Z <- matrix(data = rep(0, length = (ncol(Zinit)+nyears)*nrow(Zinit)), 
              ncol = ncol(Zinit)+nyears, nrow = nrow(Zinit), dimnames = list(dimnames(Zinit)[[1]]))
  Z[1:nrow(Zinit), 1:ncol(Zinit)] <- Zinit
  
  # initialize number of calves from Zinit
  ncalves <- length(which(rowSums(Zinit) == 1 & Zinit[,ncol(Zinit)] == 1))
  ifelse(length(ncalves) == 0, ncalves <- 0, ncalves <- ncalves)
  
  for (t in ncol(Zinit):(ncol(Zinit)+nyears-1)){ # loop from last year of Zinit
    
    nalive <- sum(Z[,t]) - ncalves # 1+ at time t
    ncalves <- 0 # reset calf counter
    ft <- f0 + (fmax-f0)*(1-(nalive/K)^z) # fec at t
    # Keep ft from becoming negative 
    # TODO: is this the correct formulation for ft?
    ifelse(ft<0, ft <- 0, ft <- ft)
    
    for (i in 1:nrow(Z)){ # loop over each individual
      
      if(Z[i, t] != 1) next else{ # if the individual isn't alive, skip
        
        if(sum(Z[i,])==51){Z[i, (t+1):ncol(Z)] <- 0} else { # if 50, doesn't survive
          
        isurvive <- rbinom(n = 1, size = 1, prob = S[sum(Z[i,])])
        ifelse(isurvive == 1, Z[i, t+1] <- 1, Z[i, (t+1):ncol(Z)] <- 0)
        
        if(isurvive == 1 & AFR < sum(Z[i,])) { # if of repro age, calf poss
          ncalves <- ncalves + rbinom(n = 1, size = 1, prob = ft)}
      
        }} # end elses
      } # end for i
    
    if (ncalves == 0) next else{
      
    # add calves to the bottom
    newrows = matrix(data = rep(0, times = ncalves*ncol(Z)), nrow = ncalves,
                     dimnames = list(stri_rand_strings(n = ncalves, length = 10)))
    newrows[,(t+1)] <- 1
    Z <- rbind(Z, newrows)} # end else
    
  } # end for t
  
  return(Z)
  
} # end function

