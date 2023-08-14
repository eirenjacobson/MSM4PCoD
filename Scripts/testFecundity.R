

K <- 100
nalive <- 0:300
S0 <- 0.8
S1 <- 0.85
S2 <- 0.95
f0 <- (1-S2) / ( (S0^AJU) * (S1^(ASA-AJU)) * (S2^(AFR-ASA)) * (1-(S2^(AMAX-AFR))) )
f0 <- 0.159
fmax <- 0.2
z <- 2.39

ft <- f0 + (fmax-f0)*(1-(nalive/K)^z) 

ft <- f0 + (fmax-f0)*(1-(nalive/K)^z)




plot(nalive/K, ft)
abline(h=0.159)
