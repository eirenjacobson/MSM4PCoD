
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(coda)

# NOTE: currently must manually create folder in Figures/ to match id

id <- "D50_LCP_bluewhale_2023-05-26"
nsim <- 10


#########################

load(paste0("./Data/MSM4PCoD_Results_", id, ".RData"))
load(paste0("./Data/MSM4PCoD_SimData_", id, ".RData"))

folder <- paste0("./Figures/", id)



#true.df <- data.frame("K1" = 100, "K2" = 50, "PCap" = 0.2, "S2" = 0.95)

# Ntot and K model results
rdf <- data.frame()
# simulated population (A + B)
sdf <- data.frame()
# region B info from simulation
bdf <- data.frame()
# lt data from simulation
ltdf <- data.frame()
# pam data from simulation
pamdf <- data.frame()
# all samples for monitored pars
pardf <- data.frame()

# create a data frame with results from all iterations
for (i in 1:nsim){
  
  r <- results[[i]]
  c1 <- data.frame(r$chain1, "Chain" = 1, "Sample" = 1:nrow(r$chain1))
  c2 <- data.frame(r$chain2, "Chain" = 2, "Sample" = 1:nrow(r$chain2))
  c3 <- data.frame(r$chain3, "Chain" = 3, "Sample" = 1:nrow(r$chain3))
  c4 <- data.frame(r$chain4, "Chain" = 4, "Sample" = 1:nrow(r$chain4))
  
  samples <- rbind.data.frame(c1, c2, c3, c4) %>% 
    select(K1, K2, S2, PCap, Chain, Sample) %>%
    mutate(Iter = i)
  
  pardf <- rbind.data.frame(pardf, samples)
  
  qdf <- summary(results[[i]])$quantiles
  ndf <- data.frame(Year = 1:100,
                    Ntot = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                               which(colnames(qdf)=="50%")],
                    LCI = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                              which(colnames(qdf)=="2.5%")],
                    UCI = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                              which(colnames(qdf)=="97.5%")],
                    K = c(rep(qdf[which(rownames(qdf) == "K1"), 
                                  which(colnames(qdf)=="50%")], 50),
                          rep(qdf[which(rownames(qdf) == "K2"), 
                                  which(colnames(qdf)=="50%")], 50)),
                    KLCI = c(rep(qdf[which(rownames(qdf)=="K1"),
                                     which(colnames(qdf)=="2.5%")], 50),
                             rep(qdf[which(rownames(qdf)=="K2"),
                                     which(colnames(qdf)=="2.5%")],50)),
                    KUCI = c(rep(qdf[which(rownames(qdf)=="K1"),
                                     which(colnames(qdf)=="97.5%")],50),
                             rep(qdf[which(rownames(qdf)=="K2"),
                                     which(colnames(qdf)=="97.5%")], 50)))
  
  rdf <- rbind.data.frame(rdf, data.frame("Iter" = i, ndf))
  
  
  sdf <- rbind.data.frame(sdf, simdata[[i]]$NSim %>%
                            pivot_wider(names_from = Region, values_from = N) %>%
                            mutate(Ntot = A + B) %>% mutate(Iter = i) %>% select(Iter, Year, Ntot))
  
  
  bdf <- rbind.data.frame(bdf, simdata[[i]]$NSim %>% filter(Region == "B") %>% 
                            mutate(Iter = i) %>% select(Iter, Year, N))
  
  CVLT <- 0.6
  CLT <- exp(1.96 * sqrt(log(1+CVLT^2)))
  ltdf <- rbind.data.frame(ltdf, simdata[[i]]$LTData %>% 
                             mutate(Iter = i) %>%
                             mutate(LCI = Nhat/CLT) %>%
                             mutate(UCI = Nhat*CLT))
  
  CVPAM <- 0.356
  CPAM <- exp(1.96 * sqrt(log(1+CVPAM^2)))
  pamdf <- rbind.data.frame(pamdf, simdata[[i]]$PAMData %>% 
                              mutate(Iter = i) %>%
                              mutate(LCI = Nhat/CVPAM) %>%
                              mutate(UCI = Nhat*CVPAM))
  
} # end for i

results.out <- list(pardf, qdf, rdf, ndf, rdf, sdf, bdf, ltdf, pamdf)

save(results.out, file = paste0("./Results/", id, ".RData"))


# now plot each iteration separately
for (i in 1:nsim){
  rtest <- filter(rdf, Iter == i)
  stest <- filter(sdf, Iter == i)
  ltest <- filter(ltdf, Iter == i)
  ptest <- filter(pamdf, Iter == i)
  btest <- filter(bdf, Iter == i)

# Plot at the level of the metapopulation (A + B)
  metaplot <- ggplot()+
                geom_ribbon(data = rtest, aes(x=Year, ymin = KLCI, ymax = KUCI), fill = "lightblue")+
                geom_line(data = rtest, aes(x=Year, y = K), color = "mediumblue", lty = "dashed")+
                geom_line(data = stest, aes(x=Year, y = Ntot, group = Iter))+
                geom_point(data = ltest, aes(x=Year, y = Nhat, group = Iter))+
                geom_errorbar(data = ltest, aes(x=Year, ymin = LCI, ymax = UCI, group = Iter))+
                geom_ribbon(data = rtest, aes(x=Year, ymin = LCI, ymax = UCI), fill = "orange2", alpha = 0.5) +
                geom_line(data = rtest, aes(x=Year, y = Ntot, group = Iter), color = "orange2") +
                ggtitle("Metapopulation") +
                theme_bw()
  
  #ggsave(metaplot, filename = paste0(folder, "/metaplot_", i, ".png"), width = 5, height = 3, units = "in")

# Plot at the level of the subpopulation (B only)

  subplot <- ggplot()+
              geom_ribbon(data = rtest, aes(x=Year, ymin = KLCI/2, ymax = KUCI/2), fill = "lightblue", alpha = 0.5)+
              geom_line(data = rtest, aes(x=Year, y = K/2), color = "mediumblue", lty = "dashed")+
              geom_line(data = btest, aes(x=Year, y = N))+
              geom_point(data = ptest, aes(x=Year, y = Nhat, group = Iter), pch = 17)+
              geom_errorbar(data = ptest, aes(x=Year, ymin = LCI, ymax = UCI, group = Iter))+
    geom_ribbon(data = rtest, aes(x=Year, ymin = LCI/2, ymax = UCI/2), fill = "orange2", alpha = 0.5) +
    geom_line(data = rtest, aes(x=Year, y = Ntot/2, group = Iter), color = "orange2", alpha = 0.5) +
              ggtitle("Subpopulation") +
              theme_bw()
  
  #ggsave(subplot, filename = paste0(folder, "/subplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  # Plot of parameters to be estimated
  
  s <- pardf %>% filter(Iter == i) %>%
    pivot_longer(cols = c(K1, K2, S2, PCap), names_to = "Parameter", values_to = "Value")
  
  parplot <- ggplot(s)+
    geom_density(aes(x=Value, color = Iter), color = "mediumpurple", fill = "mediumpurple", alpha = 0.5)+
    facet_wrap(~Parameter, scales = "free") +
    theme_bw()
  
  #ggsave(parplot, filename = paste0(folder, "/parplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  chainplot <- ggplot(s) +
    geom_path(aes(x = Sample, y= Value, color = factor(Chain)), show.legend=FALSE)+
    facet_wrap(~Parameter, scales = "free") +
    theme_bw()
  
  #ggsave(chainplot, filename = paste0(folder, "/chainplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  ggsave(ggarrange(metaplot, subplot, chainplot, parplot, nrow = 2, ncol = 2), 
         filename = paste0(folder, "/", id, "_", i, ".png"), 
         width = 12, height = 8, units = "in")

} # end for i
