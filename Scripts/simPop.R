
# Function to generate Ziphiid-ish population
# Operates deterministically by COHORTS
# Last modified by EKJ 2023-08-10

simPop <- function(S0 = 0.85, S1 = 0.9, S2 = 0.95, 
                   AFR = 10, AJU = 3, ASA = 5, AMAX = 50, 
                   fmax = 0.2, K = 1000, nyears = 1000){
    
    # S0 = calf survival
    # S1 = juvenile survival
    # S2 = subadult and adult survival
    # AFR = age at first reproduction
    # AJU = age at which become juvenile (subject to S1)
    # ASA = age at which become subadult (subject to S2)
    # fmax = maximum fecundity
    # K = 1+ carrying capacity

    z <- 2.39 # degree of compensation
    
    # construct a vector of survival probabilities
    S <- rep(NA, AMAX)
    S[1:AJU] <- S0
    S[(AJU+1):ASA] <- S1
    S[(ASA+1):AMAX] <- S2
    S[AMAX] <- 0
  
    # Calculate fec at equilibrium according to Eq 3 in Brandon et al
    # Density-dependence affects fecundity according to Pella-Tomlinson
    f0 <- (1-S2) / ( (S0^AJU) * (S1^(ASA-AJU)) * (S2^(AFR-ASA)) * (1-(S2^(AMAX-AFR))) )
    # f0 <-  0.1890236
    N <- matrix(rep(NA, nyears*(AMAX+1)), ncol=(AMAX+1)) # population matrix

    P <- rep(NA, (AMAX+1)) 
    P[1] <- 1

    for(i in 2:(AMAX+1)){
      P[i] <- prod(S[1:(i-1)])
    }
   
    R0 <- K/sum(P[2:length(P)])
    N[1, ] <- P * R0
      
    ft <-vector()
    for (t in 1:(nyears-1)){
      ft[t] <- (f0 + (fmax-f0)*(1-(sum(N[t, 2:(AMAX+1)])/K)^z)) # fec at t
      N[t+1, 1] <- sum(N[t, (AFR+1):(AMAX+1)])*ft[t] # calves
      for (a in 2:(AMAX+1)){
        N[t+1, a] <- N[t, a-1]*S[a-1] # juv and adult 
      } # end for a
    } # end for t
      
      return(list(N=N, ft=ft))
    
  } # end function

 # Ndf <- as.data.frame(simPop()$N)
 # plot(rowSums(Ndf), xlab = "Year", ylab = "Ntot")
 # plot(simPop()$ft)
 # 
# 
# names(Ndf) <- 1:AMAX
# Ndf$Year <- 1:nyears
# 
# Nlong <- Ndf %>% pivot_longer(cols = 1:AMAX, names_to = "Age")
# 
# ggplot(Nlong) +
#   geom_col(aes(x=Year, y = value, fill = Age)) +
#   geom_hline(yintercept = 1000)
# 
# plot(ft, xlab =)
