
procResults <- function(id, pars){  
  library(dplyr)
  library(tidyr)
  library(runjags)
  
  #id <- "D50_LCP_40Yrs_2023-07-23"
  #nsim <- 1
  # TODO read in excel spreadsheet to get info re nsim and CVs
  
  CVLT <- pars$lt_ecv
  CVPAM <- pars$pam_ecv
  nsim <- pars$nsim
  cyear <- pars$cyear
  nyears <- pars$nyears
  
  #########################
  
  load(paste0("./Data/Results_", id, ".RData"))
  load(paste0("./Data/SimData_", id, ".RData"))
  
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
  # all N samples from all iterations
  Ndf <- data.frame()
  
  # create a data frame with results from all iterations
  for (i in 1:nsim){
    
    r <- results[[i]]
    names(r) <- paste0("chain", 1:4)
    c1 <- data.frame(r$chain1, "Chain" = 1, "Sample" = 1:nrow(r$chain1))
    c2 <- data.frame(r$chain2, "Chain" = 2, "Sample" = 1:nrow(r$chain2))
    c3 <- data.frame(r$chain3, "Chain" = 3, "Sample" = 1:nrow(r$chain3))
    c4 <- data.frame(r$chain4, "Chain" = 4, "Sample" = 1:nrow(r$chain4))
    
    samples <- rbind.data.frame(c1, c2, c3, c4) %>% 
      select(K1, K2, S2, PCap, Chain, Sample) %>%
      mutate("Iter" = i)
    
    Nsamples <- rbind.data.frame(c1, c2, c3, c4) %>% 
      select(Chain, Sample, which(substr(names(c1), 1, 4) == "Ntot")) %>%
      rename_with(~gsub("Ntot.", "", .x, fixed = TRUE))  %>%
      pivot_longer(cols = paste0(1:nyears, "."), names_to = "Year", values_to = "Ntot") %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate("Iter" = i)
  
    Ndf <- rbind.data.frame(Ndf, Nsamples)
    
    pardf <- rbind.data.frame(pardf, samples)
    
    single.mcmc <- combine.mcmc(results[[i]])
    qdf <- summary(single.mcmc)$quantiles
    ndf <- data.frame(Year = 1:nyears,
                      Ntot = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                                 which(colnames(qdf)=="50%")],
                      LCI = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                                which(colnames(qdf)=="2.5%")],
                      UCI = qdf[which(substr(rownames(qdf), 1, 4) == "Ntot"),
                                which(colnames(qdf)=="97.5%")],
                      K = c(rep(qdf[which(rownames(qdf) == "K1"), 
                                    which(colnames(qdf)=="50%")], cyear),
                            rep(qdf[which(rownames(qdf) == "K2"), 
                                    which(colnames(qdf)=="50%")], nyears-cyear)),
                      KLCI = c(rep(qdf[which(rownames(qdf)=="K1"),
                                       which(colnames(qdf)=="2.5%")], cyear),
                               rep(qdf[which(rownames(qdf)=="K2"),
                                       which(colnames(qdf)=="2.5%")],nyears-cyear)),
                      KUCI = c(rep(qdf[which(rownames(qdf)=="K1"),
                                       which(colnames(qdf)=="97.5%")],cyear),
                               rep(qdf[which(rownames(qdf)=="K2"),
                                       which(colnames(qdf)=="97.5%")], nyears-cyear)))
    
    rdf <- rbind.data.frame(rdf, data.frame("Iter" = i, ndf))
    
    
    sdf <- rbind.data.frame(sdf, simdata[[i]]$NSim %>%
                              select(-K) %>%
                              pivot_wider(names_from = Region, values_from = N) %>%
                              mutate(Ntot = A + B) %>% mutate(Iter = i) %>% select(Iter, Year, Ntot))
    
    
    bdf <- rbind.data.frame(bdf, simdata[[i]]$NSim %>% filter(Region == "B") %>% 
                              mutate(Iter = i) %>% select(Iter, Year, N))
    
    
    CLT <- exp(1.96 * sqrt(log(1+CVLT^2)))
    ltdf <- rbind.data.frame(ltdf, simdata[[i]]$LTData %>% 
                               mutate(Iter = i) %>%
                               mutate(LCI = Nhat/CLT) %>%
                               mutate(UCI = Nhat*CLT))
    
    
    CPAM <- exp(1.96 * sqrt(log(1+CVPAM^2)))
    pamdf <- rbind.data.frame(pamdf, simdata[[i]]$PAMData %>% 
                                mutate(Iter = i) %>%
                                mutate(LCI = Nhat/CPAM) %>%
                                mutate(UCI = Nhat*CPAM))
    
  } # end for i
  
  results.out <- list(pardf = pardf, qdf = qdf, rdf = rdf, ndf = ndf, rdf = rdf, 
                      sdf = sdf, bdf = bdf, ltdf = ltdf, pamdf = pamdf, Ndf = Ndf)
  
  save(results.out, file = paste0("./Results/ProcResults_", id, ".RData"))
}
