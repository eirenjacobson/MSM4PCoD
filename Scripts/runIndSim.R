
source("./Scripts/initPop.R")
source("./Scripts/projPop.R")
source("./Scripts/redistributePopInd.R")

############################

nyears <- 100
set.c <- 1
caprecapyrs <- seq(10, 100, by = 5)

############################

set.seed(20221003)

# Set vector of carrying capacities
Ka <- rep(1000, nyears)
Kb <- c(rep(1000, nyears/2), rep(500, nyears/2))

#plot(x = NA, y = NA, xlim = c(0, 100), ylim = c(-10, 10))

# Set counter for number of animals in each area
Na_t <- rep(NA, nyears)
Nb_t <- rep(NA, nyears)

for (t in 1:nyears){
  # if it's the first year, initialize the population
  if (t == 1){
    Za_t <- projPop(Zinit = initPop(), nyears = 1)
    Zb_t <- projPop(Zinit = initPop(), nyears = 1)
  } else {
    Za_t <- projPop(Zinit = Za_tplus1, nyears = 1, K = Ka[t])
    Zb_t <- projPop(Zinit = Zb_tplus1, nyears = 1, K = Kb[t])
  }
  
  Z_new <- redistributePop(Za = Za_t, Zb = Zb_t, Ka = Ka[t], Kb = Kb[t], c = set.c)
  # These z-matrices will be the starting point for the next iteration
  Za_tplus1 <- Z_new$Za_new
  Zb_tplus1 <- Z_new$Zb_new
  
  # Fill in vectors with total number of animals in each area
  Na_t[t] <- sum(Za_tplus1[,ncol(Za_tplus1)])
  Nb_t[t] <- sum(Zb_tplus1[,ncol(Zb_tplus1)])
  
  # now simulate surveys
  if (t %in% caprecapyrs){
    
    if(t == caprecapyrs[1]){
      caprecap <- data.frame("ID" =simCapRecap(Zmat = Zb_tplus1, Pcap = 0.1),  
                             "Year" = t, "Cap" = 1)
    } else {
    caprecap <- rbind.data.frame(caprecap, data.frame("ID" =simCapRecap(Zmat = Zb_tplus1, Pcap = 0.1),  
                                                      "Year" = t, "Cap" = 1))}
  } # end caprecap
  
  if (t %in% linetransyrs){
    
    
    
  }
  
} # end for t

caphist <- caprecap %>% pivot_wider(id_cols = ID, names_from = Year, values_from = Cap)


