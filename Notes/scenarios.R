
# aim for 20 simulations

Model <- c(1, 2) # pop structure known/unknown

# add as only the version with CV of 0.6
LineTransectCV <- c(0, 0.2, 0.3, 0.6)

CapRecapEff <- c(0, 1, 2, 3)

# megan to send est of 
AcousticEff <- c(0, 1)

DeltaK <- c(0, -.25, -.5)

expand.grid(Model, LineTransectCV, CapRecapEff, AcousticEff, DeltaK)

# 256 scenarios