
procResults <- function(id, CVLT, CVPAM){  
  library(dplyr)
  library(tidyr)
  library(runjags)
  
  #id <- "NULL_LCP_Ideal1_2023-08-09"
  #nsim <- 100
  # TODO read in excel spreadsheet to get info re nsim and CVs
  
  #########################
  
  load(paste0("./Data/Results_", id, ".RData"))
  load(paste0("./Data/SimData_", id, ".RData"))
  
  # Ntot and K model results
  rdf <- data.frame()
  mdf <- data.frame()
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
  
  cdf <- data.frame()
  ncdf <- data.frame()
  
  # create a data frame with results from all iterations
  for (i in 1:nsim){
    
    r <- results[[i]]
    names(r) <- paste0("chain", 1:4)
    c1 <- data.frame(r$chain1, "Chain" = 1, "Sample" = 1:nrow(r$chain1))
    c2 <- data.frame(r$chain2, "Chain" = 2, "Sample" = 1:nrow(r$chain2))
    c3 <- data.frame(r$chain3, "Chain" = 3, "Sample" = 1:nrow(r$chain3))
    c4 <- data.frame(r$chain4, "Chain" = 4, "Sample" = 1:nrow(r$chain4))
    
    fsamples <- rbind.data.frame(c1, c2, c3, c4) %>%
      select(Chain, Sample, which(substr(names(c1), 1, 2) == "ft")) %>%
      rename_with(~gsub("ft.", "", .x, fixed = TRUE)) %>%
      pivot_longer(cols = paste0(1:(nyears-1), "."), names_to = "Year", values_to = "ft") %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate(Iter = i)
    
    nsamples <- rbind.data.frame(c1, c2, c3, c4) %>%
      select(Chain, Sample, which(substr(names(c1), 1, 2) == "ft")) %>%
      rename_with(~gsub("ft.", "", .x, fixed = TRUE)) %>%
      pivot_longer(cols = paste0(1:(nyears-1), "."), names_to = "Year", values_to = "ft") %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate(Iter = i)
    
    samples <- rbind.data.frame(c1, c2, c3, c4) %>% 
      select(K1, K2, S2, PCap, Chain, Sample) %>%
      mutate("Iter" = i)
    
    Nsamples <- rbind.data.frame(c1, c2, c3, c4) %>% 
      select(Chain, Sample, which(substr(names(c1), 1, 4) == "Ntot")) %>%
      rename_with(~gsub("Ntot.", "", .x, fixed = TRUE))  %>%
      pivot_longer(cols = paste0(1:nyears, "."), names_to = "Year", values_to = "Ntot") %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate("Iter" = i)
    
    Nncalves <- rbind.data.frame(c1, c2, c3, c4) %>% 
      select(Chain, Sample, which(substr(names(c1), 1, 9) == "noncalves")) %>%
      rename_with(~gsub("noncalves.", "", .x, fixed = TRUE)) %>%
      pivot_longer(cols = paste0(1:99, "."), names_to = "Year", values_to = "Noncalves") %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate("Iter" = i)
    
    ftN <- left_join(Nsamples, fsamples, by = c("Chain", "Sample", "Year", "Iter"))
    
    K1vals <- select(samples, K1, Chain, Sample, Iter)
    K2vals <- select(samples, K2, Chain, Sample, Iter)
    
    KvN <- left_join(Nncalves, K1vals, by = c("Chain", "Sample", "Iter")) %>%
      left_join(K2vals, by = c("Chain", "Sample", "Iter"))
    
    ggplot()+
      geom_point(data = filter(KvN, Year > 50), aes(x=Noncalves, y = K2)) +
      geom_line(aes(x=170:195, y = 170:195), col = "red")
    
    
    
    d1 <- filter(ftN, Year %in% 1:cyear) %>% left_join(K1vals, by = c("Chain", "Sample", "Iter")) %>%
      left_join(filter(Nncalves, Year %in% 1:cyear), by = c("Chain", "Sample", "Year", "Iter"))
    
    d2 <- filter(ftN, Year %in% cyear:nyears) %>% left_join(K2vals, by = c("Chain", "Sample", "Iter"))
    
    d1$NK <- d1$Noncalves/d1$K1
    d2$NK <- d2$Ntot/d2$K2
    
    d1 %>%
      group_by(Chain, Sample, Iter) %>%
    ggplot() +
      geom_line(aes(x=Year, y = NK))
    
    
    plot(nalive/K, ft, col = "red")
    
    points(d1$NK, d1$ft)
    points(d2$NK, d2$ft)
  
    Ndf <- rbind.data.frame(Ndf, Nsamples)
    
    mean.ndf <- Ndf %>% filter(Iter == i) %>% group_by(Year) %>%
      summarize(Mean = mean(Ntot), .groups = "keep")
    
    pardf <- rbind.data.frame(pardf, samples)
    
    single.mcmc <- combine.mcmc(results[[i]])
    qdf <- summary(single.mcmc)$quantiles
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
    mdf <- rbind.data.frame(mdf, cbind("Iter" = i, mean.ndf))
    
    sdf <- rbind.data.frame(sdf, simdata[[i]]$NSim %>%
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
                      sdf = sdf, bdf = bdf, ltdf = ltdf, pamdf = pamdf, Ndf = Ndf, mdf = mdf)
  
  save(results.out, file = paste0("./Results/ProcResultsMEAN_", id, ".RData"))
}
