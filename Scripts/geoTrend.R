
library(dplyr)
sdf <- results.out$sdf
ltdf <- results.out$ltdf
pamdf <- results.out$pamdf
rdf <- results.out$rdf
Ndf <- results.out$Ndf

trend <- data.frame("Iter" = rep(NA, 100),
                    "Pre" = rep(NA, 100),
                    "Post" = rep(NA, 100))
for (i in 1:100){
  for (j in unique(Ndf$Chain)){
    for (k in 1:250){
  
  predf <- filter(Ndf, Iter == i, Chain == j, Sample == k, Year <= 50)
  postdf <- filter(Ndf, Iter == i, Chain == j, Sample == k, Year > 50)

  
  
    }}}




geotrend <- function(x){
  (x[length(x)]/x[1])^(1/(length(x)-1)) 
}

geotrend(predf$Ntot)
geotrend(postdf$Ntot)

m1 <- glm(Ntot ~ Year, 
          data = predf,
          family = "quasipoisson")

m2 <- lm(log(Ntot) ~ Year, data = filter(Ndf, Year %in% 1:50, Iter == i, C hain == j, Sample == k))

m3 <- 

pred <- data.frame("Year" = 1:50, "Pred" = NA)
pred$Pred <- log(predict(m1, newdata = pred, type = "response"))
pred$GeoPred <- log(predf$Ntot[1]*geotrend(predf$Ntot)^(pred$Year))

 



