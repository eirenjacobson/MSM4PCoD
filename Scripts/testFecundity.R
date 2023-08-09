

K <- 200
nalive <- 100:300
S0 <- 0.8
S1 <- 0.85
S2 <- 0.95
f0 <- (1-S2)/(S0^(AJU-1)*S1^(ASA-1)*S2^(AFR-ASA-2))*(1-S2^(AMAX-AFR-2))
fmax <- 0.2
z <- 2.39

ft <- f0 + (fmax-f0)*(1-(nalive/K)^z) 

ft <- f0 + (fmax-f0)*(1-(nalive/K)^z)

          