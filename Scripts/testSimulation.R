
  ####
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nimble)
  library(nimbleEcology)
  ####
  
  source("./Scripts/initPop.R")
  source("./Scripts/projPop.R")
  source("./Scripts/redistributePopInd.R")
  source("./Scripts/simCapRecap.R")
  source("./Scripts/simSurvey.R")
  
  ############################

  # years to run simulation
  nyears <- 100
  # connectivity parameter
  cval <- 1

  Ka_1 <- 100
  Kb_1 <- 100
  Ka_2 <- 100
  Kb_2 <- 100

  
  ############################
  
  # Set vector of carrying capacities
  Ka <- c(rep(Ka_1, nyears/2), rep(Ka_2, nyears/2))
  Kb <- c(rep(Kb_1, nyears/2), rep(Kb_2, nyears/2))
  
  # Set counter for number of animals in each area
  Na_t <- rep(NA, nyears)
  Nb_t <- rep(NA, nyears)
  
  for (t in 1:nyears){
    # if it's the first year, initialize the population
    if (t == 1){
      Za_t <- projPop(Zinit = initPop(K=Ka_1), nyears = 1)
      Zb_t <- projPop(Zinit = initPop(K=Kb_1), nyears = 1)
    } else {
      Za_t <- projPop(Zinit = Za_tplus1, nyears = 1, K = Ka[t])
      Zb_t <- projPop(Zinit = Zb_tplus1, nyears = 1, K = Kb[t])
    }
    
    Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t], c = cval)
    # These z-matrices will be the starting point for the next iteration
    Za_tplus1 <- Z_new$Za_new
    Zb_tplus1 <- Z_new$Zb_new
    
    # Fill in vectors with total number of animals in each area
    Na_t[t] <- sum(Za_tplus1[,ncol(Za_tplus1)])
    Nb_t[t] <- sum(Zb_tplus1[,ncol(Zb_tplus1)])
    
    
  } # end for t
  
  N <- data.frame("Year" = rep(1:nyears, 2),
                  "N" = c(Na_t, Nb_t), 
                  "Region" = c(rep("A", length(Na_t)), rep("B", length(Nb_t))),
                  "K" = c(Ka, Kb))
  
  ggplot(N) +
    geom_line(aes(x=Year, y = N, color = Region))
  
  
  N %>% pivot_wider(id_cols = c(Year), names_from = Region, values_from = N) %>%
    mutate(Ntot = A + B) %>%
    ggplot()+
    geom_line(aes(x=Year, y = Ntot))
  