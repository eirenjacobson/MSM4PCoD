
source("./Scripts/NimbleIPM_NoPopModel.R")

set.seed(20230925)
trueK <- 500

LT.est <- rnorm(100, mean = trueK, sd = 50)
PAM.data <- rnorm(100, mean = trueK, sd = 50)

nimbleData <- list(ltestN = LT.data, pamestN = PAM.data)

nimbleConstants <- list(K_lower = 400,
                        K_upper = 600)

nimbleInits <- list(K1_scalar = 1)

nimbleParams <- list("K1")

model <- nimbleModel(code = ipm,
                     constants = nimbleConstants,
                     data = nimbleData,
                     inits = nimbleInits,
                     check = FALSE)

# Run the model
nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        constants = nimbleConstants, data = nimbleData,
                        thin = 10, niter = 10000, nburnin = 9000, nchains = 4,
                        samplesAsCodaMCMC = TRUE)

summary(nimbleOut)
