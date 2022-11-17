
simLineTransect <- function(N, L, A, ESHW = 2, Tdist = 4, g0 = 0.55){
  
  # Inputs
  # N: total (true) number of animals in the study area
  # L: trackline distance covered (km)
  # A: study area size (km2)
  # ESHW: effective strip half-width (km) default 2km
  # g0: trackline detection probability, default 0.55
  
  EES <- L*ESHW*2 # effective area searched
  pEES <- EES/A # prob of being in the EES

  Nd <- rbinom(n = N, size = 1, prob = pEES*g0)
 
  Nhat <- sum(Nd)*A/EES*(1/g0)
  
  SD <- sqrt(N*pEES*g0*(1-(pEES*g0)))
 
  LTEst <- list(Nhat = Nhat, SD = SD)
  
  return(LTEst)
  
}



