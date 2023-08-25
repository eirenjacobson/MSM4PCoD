
library(readxl)
library(dplyr)
library(tidyr)
filepath <- "./Data/MSM4PCoD_SimulationParameters_V3.xlsx"
nscenarios <- 11
params <- read_excel(filepath, 
                     col_types = c("text", rep("numeric", nscenarios)), 
                     na = "NA")
scenarios <- c("NULL_Ideal",    "D50B_Real",   "D75B_Real",   "NULL_Real")  
d <- "2023-08-21"

for (i in 1:length(scenarios)){


pars <- params[,c(1,which(names(params) == scenarios[i]))] %>% 
  pivot_wider(names_from = "...1", values_from = scenarios[i])

id <- paste0(scenarios[i], "_", d)

source("./Scripts/procResults.R")
procResults(id, pars)

source("./Scripts/calcPower_v2.R")
calcPower(id)
}