
# check lognormal transformation
# example target N and CV
  N <- 200
  CV <- 0.3

# the target mean and sd on the real scale
  mu_X <- N # mean
  sigma_X <- CV*N # standard deviation
  sigma2_X <- (sigma_X)^2 # variance

# transform to mean and sd on the log scale
  mu <- log(mu_X^2/sqrt(mu_X^2 + sigma2_X))
  sigma2 <- log(1+(sigma2_X/mu_X^2)) #variance
  sigma <- sqrt(sigma2)

# check that coverage matches expected CIs
  C <- exp(1.96 * sqrt(log(1+CV^2)))
  LCI <- N/C
  UCI <- N*C
  samples <- rlnorm(10000, meanlog = mu, sdlog = sigma)
  length(samples[which(samples>LCI & samples<UCI)])/length(samples)
  mean(samples)
  median(samples)
 
 # draw some random "observations"
 set.seed(20230928)
 
  NEst <- rlnorm(1000, meanlog = mu, sdlog = sigma)  
 
 # now test the likelihood of those observations given different
 # model estimates of N (these are the i values in the loop)
 # convert each observation to the log scale
  mu_NEst <- NEst
  sigma_NEst <- CV*N
  
  NEst_log <- log(mu_NEst^2/sqrt(mu_NEst^2 + sigma_NEst^2))
  NEst_sigma_log <- sqrt(log(1+(sigma_NEst^2/mu_NEst^2)))
  
  i_log <- log(mu_NEst^2/sqrt(mu_NEst^2 + sigma_NEst^2))
  i_sigma <- sqrt(log(1+(sigma_NEst^2/i^2)))
 
  LL <- vector()
  
  for (i in 150:250){
    LL[i] <- sum(dlnorm(NEst, log(i), sdlog = NEst_sigma_log))
  }
  
  plot(150:250, LL[150:250])
  points(which.max(LL), LL[which.max(LL)], col = "red", pch = 16)
  abline(v = 200)
 
  