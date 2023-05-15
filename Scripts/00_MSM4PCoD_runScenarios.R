


##############
# path to file containing simulation parameters
filepath <- "./Data/MSM4PCoD_SimulationParameters.xlsx"
# number of scenarios IN THE FILE
nscenarios <- 2
# names of the scenarios you want to run
scenarios <- c("A", "B")
# set seed for reproducibility
set.seed(20230504)
#TODO: add options to runIndSim and nimbleIPM to include/not include survey types
#TODO: change read_excel to only read columns of interest (pars + scenarios)
##############

start <- Sys.time()

library(readxl)
library(dplyr)
library(tidyr)

source("./Scripts/runIndSim.R")

params <- read_excel(filepath, 
                     col_types = c("text", rep("numeric", nscenarios)), 
                     na = "NA")
#TODO: structure results so it can accommodate multiple scenarios/iterations
results <- list()
for (i in 1:length(scenarios)){
  
  pars <- params[,c(1,which(names(params) == scenarios[i]))] %>% 
    pivot_wider(names_from = "...1", values_from = scenarios[i])
  
  if(pars$linetrans == TRUE){
    linetransyrs <- seq(from = pars$lt_start, 
                        to = pars$lt_end, 
                        by = pars$lt_int)}
  if(pars$caprecap == TRUE){
    caprecapyrs <- seq(from = pars$cr_start, 
                       to = pars$cr_end, 
                       by = pars$cr_int)}
  
  for (j in 1:pars$nsim){
    print(paste("Beginning iteration", j, "of scenario", scenarios[i]))
    results[i][[j]] <- runIndSim(nyears = pars$nyears, 
                              cval = pars$cval,
                              Ka_1 = pars$Ka_1,
                              Ka_2 = pars$Ka_2,
                              Kb_1 = pars$Kb_1, 
                              Kb_2 = pars$Kb_2,
                              linetransyrs = linetransyrs, 
                              lt_ecv = pars$lt_ecv, 
                              caprecapyrs = caprecapyrs, 
                              pcap = pars$pcap)
    
  } # end for j
} # end for i

end <- Sys.time()
dur <- end-start 

save(results, file = "./Data/results_scenarioB.RData")
