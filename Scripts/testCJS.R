library(nimble)
library(nimbleEcology)

load("./Data/MSM4PCoD_SimData_D50_LCP_2023-06-27.RData")

source("./Scripts/cjs.R")

#for (i in 1:100){

i <- 2
sd <- simdata[[i]]

cr <- sd$CRData
crstart <- min(cr$Year)
crend <- max(cr$Year)
ncryears <- length(crstart:crend)

caphist <- complete(cr, ID, Year = crstart:crend, fill = list(Cap = 0)) %>% 
  pivot_wider(id_cols = ID, names_from = Year, values_from = Cap)

Y <- as.matrix(caphist[,2:ncol(caphist)])
# take out individuals only seen on last occasion
lastind <- which(rowSums(Y) == 1 & Y[,ncol(Y)] == 1)
if(length(lastind) == 0){Y <- Y} else {Y <- Y[-lastind,]}

Find <- rep(0, nrow(Y))
for (i in 1:nrow(Y)){
  Find[i] <- min(which(Y[i,] == 1))
}

nimbleData <- list(Y = Y)
nimbleInits <- list("S2" = 0.95, "PCap" = 0.2)
nimbleConstants <- list(ncryears = ncryears, Nind = nrow(Y), Find = Find)
nimbleParams <- list("S2", "PCap")

model <- nimbleModel(code = cjs,
                     data = nimbleData,
                     constants = nimbleConstants,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        constants = nimbleConstants, data = nimbleData,
                        thin = 10, niter = 10000, nburnin = 9000, nchains = 4,
                        samplesAsCodaMCMC = TRUE)

summary(nimbleOut)

