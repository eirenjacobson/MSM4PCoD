

##############
# path to file containing simulation parameters
filepath <- "./Data/MSM4PCoD_SimulationParameters.xlsx"
# number of scenarios IN THE FILE
nscenarios <- 7
# names of the scenarios you want to run
scenarios <- c("NULL_LCP", "D50_LCP")
# set seed for reproducibility
set.seed(20230504)
##############

start <- Sys.time()

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)

source("./Scripts/runIndSim.R")
source("./Scripts/runNimble.R")

params <- read_excel(filepath, 
                     col_types = c("text", rep("numeric", nscenarios)), 
                     na = "NA")
simdata <- list()
results <- list()
for (i in 1:length(scenarios)){
  
  pars <- params[,c(1,which(names(params) == scenarios[i]))] %>% 
    pivot_wider(names_from = "...1", values_from = scenarios[i])
  
  if(pars$linetrans == TRUE){
    linetransyrs <- seq(from = pars$lt_start, 
                        to = pars$lt_end, 
                        by = pars$lt_int)} else {linetransyrs <- NA}
  if(pars$caprecap == TRUE){
    caprecapyrs <- seq(from = pars$cr_start, 
                       to = pars$cr_end, 
                       by = pars$cr_int)} else {capreapyrs <- NA}
  if(pars$pam == TRUE){
    pamyrs <- seq(from = pars$pam_start,
                  to = pars$pam_end,
                  by = pars$pam_int)} else {pamyrs <- NA}
  
  for (j in 1:pars$nsim){
    print(paste("Beginning iteration", j, "of scenario", scenarios[i]))
    
    simdata[[j]] <- runIndSim(nyears = pars$nyears, 
                        cval = pars$cval,
                        Ka_1 = pars$Ka_1,
                        Ka_2 = pars$Ka_2,
                        Kb_1 = pars$Kb_1, 
                        Kb_2 = pars$Kb_2,
                        linetrans = pars$linetrans,
                        linetransyrs = linetransyrs, 
                        lt_ecv = pars$lt_ecv, 
                        caprecap = pars$caprecap,
                        caprecapyrs = caprecapyrs, 
                        pcap = pars$pcap,
                        pam = pars$pam,
                        pamyrs = pamyrs,
                        pam_ecv = pars$pam_ecv)
                        
    results[[j]] <- runNimble(simdata = simdata[[j]],
                              linetrans = pars$linetrans,
                              caprecap = pars$caprecap,
                              pam = pars$pam,
                              nchains = pars$nchains,
                              thin = pars$thin, 
                              niter = pars$niter, 
                              nburnin = pars$nburnin)
    
  } # end for j
  
  save(results, file = paste0("./Data/MSM4PCoD_Results_Scenario", 
                              scenarios[i], "_", date(now()), ".RData"))
  
  
} # end for i

end <- Sys.time()
dur <- end-start 

