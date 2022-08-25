
# Cohort-based simulation

source("./Scripts/simPop.R")
source("./Scripts/redistributePop.R")

library(ggplot2)

library(tidyverse)

# run population simulation for 100 years to get stable pop 
Ndefault <- simPop(K = 1000, new = TRUE, nyears = 100)

# set up two neighboring populations, Na and Nb

Na <- Ndefault[100,]
Nb <- Ndefault[100,]

Ka <- rep(1000, 100)
Kb <- c(rep(1000, 50), rep(1000, 50))

Na_proj <- matrix(nrow = 100, ncol = 50)
Na_proj[1,] <- Na

Nb_proj <- matrix(nrow = 100, ncol = 50)
Nb_proj[1,] <- Nb

for (t in 1:99){
  
  # redistribute according to current population size and perceived K
  N_new <- redistributePop(Na = Na_proj[t,], Nb = Nb_proj[t,], Ka = 1000, Kb = 1000, c = 0.5)
  
  # now project forward according to actual K
  Na_proj[t+1,] <- simPop(new = FALSE, K = Ka[t], N1 = N_new$Na_new, nyears = 2)[2,]
  Nb_proj[t+1,] <- simPop(new = FALSE, K = Kb[t], N1 = N_new$Nb_new, nyears = 2)[2,]
  
}

Na_df <- as.data.frame(Na_proj)
colnames(Na_df) <- 1:50
Na_df$Year <- 1:nrow(Na_df)
Na_long <- pivot_longer(Na_df, cols = 1:50, names_to = "Age")
Na_long$Region <- "A"


Nb_df <- as.data.frame(Nb_proj)
colnames(Nb_df) <- 1:50
Nb_df$Year <- 1:nrow(Nb_df)
Nb_long <- pivot_longer(Nb_df, cols = 1:50, names_to="Age")
Nb_long$Region <- "B"

Nmat <- rbind.data.frame(Na_long, Nb_long)
Nmat$Age <- as.numeric(Nmat$Age)

ggplot(Nmat) +
  geom_point(aes(x=Year, y = Age, size = value)) +
  facet_wrap(~Region)

Nmat %>% 
  group_by(Year, Region) %>%
  summarize(Ntot = sum(value)) %>%
  ggplot() +
  geom_point(aes(x=Year, y = Ntot, color = Region)) +
  facet_wrap(~Region)


ggplot(Nmat) +
  geom_point(aes())

