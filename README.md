# MSM4PCoD

2023-08-25

- The scripts beginning 00_MSM4PCoD_runScenarios will run all code (simulation + model fitting)
- Currently set up with normal observation processes for line-transect and PAM data in both the simulation and model
- After these are run, plotScenarios.R will generate plots of the iterations of the scenarios. Note that a directory matching the scenario ID must be manually created in the Figures folder. 
- For quick testing, try changing the number of simulations (nsim in the xls file) to e.g. 10 and comment out (don’t run) the power calculation part of the routine (lines 115-116 in 00_MSM4PCoD_runScenarios)
- Note that because of file sizes/space constraints, not all results/figures have been added to this repo.

Required packages:

parallel  
nimble  
nimbleEcology  
dplyr  
tidyr  
lubridate  
readxl  
ggplot2  
foreach  
doParallel  
coda  
ggpubr  
runjags  
