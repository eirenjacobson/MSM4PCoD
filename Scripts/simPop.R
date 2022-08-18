
# Function to generate Ziphiid-ish population
# Last modified by EKJ 2022-08-17

simPop <- function(S0 = 0.8, S1 = 0.85, S2 = 0.95, 
                   AFR = 10, AJU = 3, ASA = 5, AMAX = 50,
                   fmax = 0.2, K = 1000, nyears = 100,
                   new = TRUE, N1 = NA){
    
    # S0 = calf survival
    # S1 = juvenile survival
    # S2 = subadult and adult survival
    # AFR = age at first reproduction
    # AJU = age at which become juvenile (subject to S1)
    # ASA = age at which become subadult (subject to S2)
    # fmax = maximum fecundity
    # K = carrying capacity
    # new = whether or not to generate a new population
    # if new = FALSE, N1 is the row of 

    z <- 2.39 # degree of compensation
    
    # construct a vector of survival probabilities
    S <- rep(NA, AMAX-1)
    S[1:(AJU-1)] <- S0
    S[AJU:(ASA-1)] <- S1
    S[ASA:AMAX] <- S2
  
    # Calculate fec at equilibrium according to Eq 3 in Brandon et al
    # Density-dependence affects fecundity according to Pella-Tomlinson
    f0 <- (1-S2)/(S0^(AJU-1)*S1^(ASA-1)*S2^(AFR-ASA-2))*(1-S2^(AMAX-AFR-2))
    
    if (new == TRUE){
      
      N <- matrix(rep(NA, nyears*AMAX), ncol=AMAX) # population matrix

      P <- rep(NA, x) 
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
      N[1, ] <- P * R0
    
    } # end if new
    
    if (new == FALSE) {
      
      N <- matrix(rep(NA, (nyears+1)*AMAX), ncol=AMAX)
      N[1,] <- N1
      
    } # end if else
    
    for (t in 1:(nyears-1)){
      ft <- (f0 + (fmax-f0)*(1-(sum(N[t, 2:AMAX])/K)^z)) # fec at t
      print(ft)
      N[t+1, 1] <- sum(N[t, AFR:AMAX])*ft # calves
      for (a in 2:AMAX){
        N[t+1, a] <- N[t, a-1]*S[a-1] # juv and adult 
      } # end for a
    } # end for t
      
      return(N)
    
  } # end function

