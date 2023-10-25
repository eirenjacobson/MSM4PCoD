# MSM4PCoD

2023-10-25

- The scripts beginning 00_MSM4PCoD_runScenarios will run all code (simulation + model fitting)
- Currently set up with normal observation processes for line-transect and PAM data in both the simulation and model
- After these are run, plotScenarios.R will generate plots of the iterations of the scenarios (currently set to plot the first 10 iterations).
- plotPower.R will plot the power resulting from a single scenario run
- combinePower.R will plot the power and generate a table of reference power values for 2+ combined scenario runs (i.e., combining the NULL and D50AB runs for a given scenario will cover a wider range of declines than either of those on their own).
- See MSM4PCoD_Task3C_Results_2023-10-24.Rmd for an example workflow to generate power plots (once scenarios have been run).
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
