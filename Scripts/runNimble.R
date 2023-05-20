
runNimble <- function(simdata, 
                      linetrans, caprecap, pam,
                      nchains, niter, thin, nburnin){
  
  source("./Scripts/nimbleIPM.R")
  source("./Scripts/simPop.R")
  
  # Process caphist for nimbleIPM
  if(caprecap == TRUE){
    cr <- simdata$CRData
    crstart <- min(cr$Year)
    crend <- max(cr$Year)
    ncryears <- length(crstart:crend)
  
    caphist <- complete(cr, ID, Year = crstart:crend, fill = list(Cap = 0)) %>% 
      pivot_wider(id_cols = ID, names_from = Year, values_from = Cap)
  
    Y <- as.matrix(caphist[,2:ncol(caphist)])
    # take out individuals only seen on last occasion
    Y <- Y[-which(rowSums(Y) == 1 & Y[,ncol(Y)] == 1),]
  
    Find <- rep(0, nrow(Y))
    for (i in 1:nrow(Y)){
      Find[i] <- min(which(Y[i,] == 1))
    }
  } else {Y <- NULL; Find <- NULL; ncryears <- NULL} # end if caprecap
  
  Ndefault <- simPop(K = round(max(Nhat$Nhat)), new = TRUE, nyears = nyears)
  
  # stopping here -- need to figure out how to get info out of list if NULL
  # or make different versions of data list depending on which data are provided?
  nimbleData <- list(ltestN = simdata$LTData$Nhat, ltestSD = simdata$LTData$SD, 
                     Y = Y, 
                     pamestN = simdata$PAMData$Nhat, pamestSD = simdata$PAMData$SD)

    nimbleConstants <- list(cyear = 50, nyears = nyears, S0 = 0.8, S1 = 0.85, AFR = 10,
                          AJU = 3, ASA = 5, AMAX = 50, fmax = 0.25, nyears = 100, 
                          z = 2.39, ncryears = ncryears, Nind = nrow(Y),Find = Find, 
                          ltyears = simdata$LTData$Year, nltyears = length(simdata$LTData$Year),
                          pamyears = simdata$PAMData$Year, npamyears = length(simdata$PAMData$Year),
                          K1_upper = round(max(c(simdata$LTData$Nhat, simdata$PAMData$Nhat))*2), 
                          K2_upper = round(max(c(simdata$LTData$Nhat, simdata$PAMData$Nhat))*2))
  
  nimbleInits <- list(S2 = 0.95, 
                      PCap = 0.2,
                      K1 = round(max(Nhat$Nhat)), 
                      K2 = round(max(Nhat$Nhat)), 
                      N = round(Ndefault), Ntot = rowSums(Ndefault))
  
  nimbleParams <- list("S2", "K1", "K2", "Ntot", "ft", "f0")
  
  model <- nimbleModel(code = ipm,
                       constants = nimbleConstants,
                       data = nimbleData,
                       inits = nimbleInits,
                       check = FALSE)
  
  # Run the model
  nimbleOut <- nimbleMCMC(model, 
                          monitors = nimbleParams, 
                          thin = thin, niter = niter, nburnin = nburnin, nchains = nchains,
                          samplesAsCodaMCMC = TRUE)
  
  
  
  return(nimbleOut)
  
  
}