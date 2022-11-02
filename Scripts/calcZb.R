
library(readxl)
vals <- read_excel("./Data/Visual_power_parameters.xlsx")

calcZb <- function(r, n, cv, za){
  
  sqrt((r^2 * n^3)/12*(cv^2)) - za
  
}

Zb <- calcZb(r = vals$`Rate of change`, n = vals$`N value`,
                cv = vals$`Population cv`, za = 1.645)
