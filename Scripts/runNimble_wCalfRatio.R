
runNimble <- function(pars, simdata, 
                      nyears, linetrans, caprecap, pam,
                      nchains, niter, thin, nburnin, X){
  
  library(nimble)
  library(nimbleEcology)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  
  source("./Scripts/nimbleIPM_wCalfRatio.R")
  source("./Scripts/simPop.R")
  
  # Process caphist for nimbleIPM
  if(caprecap == TRUE){
    cr <- simdata$CRData
    crstart <- min(cr$Year)
    crend <- max(cr$Year)
    ncryears <- length(unique(cr$Year))
    cryears <- unique(cr$Year)
    
    caphist <- complete(cr, ID, Year = crstart:crend, fill = list(Cap = 0)) %>% 
      pivot_wider(id_cols = ID, names_from = Year, values_from = Cap)
    nadults <- colSums(caphist[,2:ncol(caphist)])
    Y <- as.matrix(caphist[,2:ncol(caphist)])
 #   # take out individuals only seen on last occasion
    lastind <- which(rowSums(Y) == 1 & Y[,ncol(Y)] == 1)
    if(length(lastind) == 0){Y <- Y} else {Y <- Y[-lastind,]}
    
    Find <- rep(0, nrow(Y))
    for (i in 1:nrow(Y)){
      Find[i] <- min(which(Y[i,] == 1))
    }
  } else {Y <- NULL; Find <- NULL; ncryears <- NULL} # end if caprecap
  
  #Kdefault <- max(simdata$LTData$Nhat, simdata$PAMData$Nhat*2)
  Ndefault <- simPop(K = 200, nyears = pars$nyears)  

  nimbleData <- list(ltestN = simdata$LTData$Nhat, ltestSDlog = simdata$LTData$SDLog, 
                     Y = Y, NCalves = simdata$CalfData$ObsCalves, NNonCalves = nadults,
                     pamestN = simdata$PAMData$Nhat, pamestSDlog = simdata$PAMData$SDLog)

    nimbleConstants <- list(cyear = pars$cyear, nyears = pars$nyears, S0 = 0.85, S1 = 0.9, AFR = 10,
                          AJU = 3, ASA = 5, AMAX = 50, fmax = 0.2,  
                          z = 2.39, ncryears = ncryears, cryears = cryears, Nind = nrow(Y), Find = Find, 
                          ltyears = simdata$LTData$Year, nltyears = length(simdata$LTData$Year),
                          pamyears = simdata$PAMData$Year, npamyears = length(simdata$PAMData$Year),
                          K_lower = 75,#round(min(c(simdata$LTData$Nhat, simdata$PAMData$Nhat*2))), 
                          K_upper = 225)#round(max(c(simdata$LTData$Nhat, simdata$PAMData$Nhat*2))))
  
  nimbleInits <- list(S2 = 0.95, 
                      PCap = 0.4,
                      K1_scalar = 1, 
                      K2_scalar = 1, 
                      N = round(Ndefault$N), Ntot = rowSums(Ndefault$N))
  
  nimbleParams <- list("S2", "K1", "K2", "PCap", "Ntot", "ft", "f0", "calves", "noncalves")
  
  model <- nimbleModel(code = ipm,
                       constants = nimbleConstants,
                       data = nimbleData,
                       inits = nimbleInits,
                       check = TRUE)

  # Run the model
  nimbleOut <- nimbleMCMC(model, 
                          monitors = nimbleParams, 
                          constants = nimbleConstants, data = nimbleData,
                          thin = thin, niter = niter, nburnin = nburnin, nchains = nchains,
                          samplesAsCodaMCMC = TRUE)
  
  return(nimbleOut)
  
  
}