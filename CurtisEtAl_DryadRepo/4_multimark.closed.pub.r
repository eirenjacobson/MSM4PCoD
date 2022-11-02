# multimark Closed

# Based on sensitivity simulations, should use no-precaps data for closed-population estimates. 
# Results for this data set support that: when precaps included, as add more years,  get lower estimates. 

library(multimark)
library(coda)
library(MCMCvis)
library(tidyverse)
library(mcmcplots)

# data setup
CH = data.matrix(read.csv("ch.yy.multi.npc.20190315.csv",colClasses="numeric"))
known = read.csv("multi.known.20190315.csv")
occ.cov = read.csv("occasions.cov.yy.20190315.csv")
eff = data.frame(eff=occ.cov$eff.bssn)
setup <- processdata(CH, data.type="sometimes", covs=eff, known=known$kn)

i <- 9:11
subCH <- CH[,i]
known.sub <- known[rowSums(subCH)>0,]
subCH <- subCH[rowSums(subCH)>0,]   # omit all-zero rows 
setup.sub <- processdata(subCH, data.type="sometimes", covs=data.frame(eff=scale(eff[i,,drop=F])), known=known.sub$kn)

# models
### RESULTS FOR P-C MODELS SHOW TRAP-SHYNESS, RATHER THAN INCREASE IN C OVER P, WHICH WOULD MAKE UNDERESTIMATION WORSE -> DROP C
# effiRE for three years
multi.effiRE.sub = multimarkClosed(mms=setup.sub,data.type="sometimes",mod.p=~1+eff+h,mod.delta=~1,nchains=3,burnin=250000,iter=1000000,thin=100,parms=c("pbeta","N","sigma2_zp","delta","alpha","psi"))
multi.effiRE.sub.p = getprobsClosed(multi.effiRE.sub)
gelman.diag(multi.effiRE.sub$mcmc)
heidel.diag(multi.effiRE.sub$mcmc)
rmeanplot(multi.effiRE.sub$mcmc)
# effiRE for all years
multi.effiRE = multimarkClosed(mms=setup,data.type="sometimes",mod.p=~1+eff+h,mod.delta=~1,nchains=3,burnin=50000,iter=200000,thin=50,parms=c("pbeta","N","sigma2_zp","delta","alpha","psi"))
multi.effiRE.p = getprobsClosed(multi.effiRE)
gelman.diag(multi.effiRE$mcmc)
rmeanplot(multi.effiRE$mcmc)

# save models
save(multi.effiRE.sub, multi.effiRE.sub.p, multi.effiRE, multi.effiRE.p, file="multimark.closed.20190315.rdata")
