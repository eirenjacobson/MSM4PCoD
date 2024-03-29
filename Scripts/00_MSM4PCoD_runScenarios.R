
library(parallel)
ncores <- 4
this_cluster <- makeCluster(ncores)

##############
# path to file containing simulation parameters
filepath <- "./Data/MSM4PCoD_SimulationParameters_V3.xlsx"
# number of scenarios IN THE FILE
nscenarios <- 21
# names of the scenarios you want to run
scenarios <- "NULL_Ideal"
# whether to simulate new data
simnewdata <- TRUE
# if simnewdata = FALSE, specify where to find datasets (for each scenario)
#simdataloc <- c("./Data/MSM4PCoD_SimData_NULL_LCP_bluewhale_2023-05-26.RData", 
#                "./Data/MSM4PCoD_SimData_D50_LCP_bluewhale_2023-05-26.RData")

# set seed for reproducibility
set.seed(20230504)
##############

start <- Sys.time()

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(parallel)

source("./Scripts/runIndSim.R")
source("./Scripts/runNimble.R")

params <- read_excel(filepath, 
                     col_types = c("text", rep("numeric", nscenarios)), 
                     na = "NA")
results <- list()
simdata <- list()

stamp <- date(now())

for (i in 1:length(scenarios)){
  
  if(simnewdata == TRUE){simdata <- list()}else{load(simdataloc[i])}
  
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
    set.seed(paste0("202305", j))
    
    if(simnewdata == TRUE){
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
                        pam_ecv = pars$pam_ecv, 
                        pars = pars)}

    results[[j]] <- parLapply(fun=runNimble, X=1:ncores, cl=this_cluster,
                              pars = pars,
                              simdata = simdata[[j]],
                              linetrans = pars$linetrans,
                              caprecap = pars$caprecap,
                              pam = pars$pam,
                              nyears = pars$nyears,
                              nchains = 1,
                              thin = pars$thin,
                              niter = pars$niter,
                              nburnin = pars$nburnin)
    
  } # end for j
  
  id <- paste0(scenarios[i], "_", stamp)

  save(results, file = paste0("./Data/Results_", id, ".RData"))
  
  save(simdata, file = paste0("./Data/SimData_", id, ".RData"))


  source("./Scripts/procResults.R")
  procResults(id, pars)

  source("./Scripts/calcPower_v2.R")
  calcPower(id, ncores)

} # end for i

stopCluster(this_cluster)

end <- Sys.time()
dur <- end-start 
