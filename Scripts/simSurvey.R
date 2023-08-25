
# generic function to simulate data collection given true N and expected CV
# returns a single survey estimate NEst drawn from a lognormal distribution

simSurvey <- function(N, CV){
  
  # N is (true/simulated) population size, assumed to be estimated without bias
  # CV is the expected CV of a survey
  
  #N <- 100
  #CV <- 0.2
  
  # the target mean and sd on the real scale
  mu_X <- N
  sigma_X <- CV*N
  
  # transform to mean and sd on the log scale
  mu <- log(mu_X^2/sqrt(mu_X^2 + sigma_X^2))
  sigma <- sqrt(log(1+(sigma_X^2/mu_X^2)))
  
  # check that coverage matches expected CIs
  # C <- exp(1.96 * sqrt(log(1+CV^2)))
  # LCI <- N/C
  # UCI <- N*C
  # samples <- rlnorm(10000, meanlog = mu, sdlog = sigma)
  # length(samples[which(samples>LCI & samples<UCI)])/length(samples)
  
  #NEst <- rlnorm(1, meanlog = mu, sdlog = sigma)
  NEst <- rnorm(1, mean = mu_X, sd = sigma_X)
  
  return(NEst)
  
}