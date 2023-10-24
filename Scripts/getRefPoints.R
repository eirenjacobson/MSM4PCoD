
# Helper functions for calculating reference points
# reorganizing N_t = N_0e^(r*t)
# N_t/N_0 = decline
# r = rate
# t = time (years)

getRate <- function(decline, time){
  
  time <- 10
  decline <- 0.5
  
  r <- log(decline)/time
  
  return(r)
  
}

getTime <- function(decline, rate){
  
  
}

getDecline <- function(rate, time){
  
  
}

r <- -0.005
time <- 10
exp(r*time)

r <- 0.005
decline <- 0.5

log(decline)*r