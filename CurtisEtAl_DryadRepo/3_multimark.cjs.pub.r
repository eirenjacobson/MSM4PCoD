# multimark CJS

library(multimark)
library(coda)
library(tidyverse)
library(mcmcplots)

# common data import
known = read.csv("multi.known.20190315.csv")
occ.cov = read.csv("occasions.cov.yy.20190315.csv")
eff = data.frame(eff02=scale(occ.cov$eff.02), eff.bssn=occ.cov$eff.bssn)

# fit when exclude "precaptures"
## data setup
CH = data.matrix(read.csv("ch.yy.multi.npc.20190315.csv",colClasses="numeric"))
setup = processdata(CH,data.type="sometimes",covs=eff,known=known$kn)
## fit models with all parameters kept
multi.dot = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1,mod.delta=~1,nchains=3,burnin=5000,iter=150000, thin=50,parms="all")
multi.iRE = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+h,mod.delta=~1,nchains=3,burnin=5000,iter=150000, thin=50,parms="all")
multi.effiRE = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff02+h,mod.delta=~1,nchains=3,burnin=5000,iter=150000, thin=50,parms="all")
multi.effbssniRE = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff.bssn+h,mod.delta=~1,nchains=3,burnin=5000,iter=150000, thin=50,parms="all")
## Compare: pdot, pire, peffire, pbssnire
multimodel.mm4 = multimodelCJS(list(multi.dot, multi.iRE, multi.effiRE, multi.effbssniRE),monparms=c("phi","p"))
gelman.diag(multimodel.mm4$rjmcmc, multivariate=F)
multimodel.mm4$rjmcmc = multimodel.mm4$rjmcmc[,c(1:10,110:111)]
## get p and phi
multi.dot.p = getprobsCJS(multi.dot)   
multi.dot.p = multi.dot.p[,c(1,110)]
multi.iRE.p = getprobsCJS(multi.iRE)
multi.iRE.p = multi.iRE.p[,c(1,110)]
multi.effiRE.p = getprobsCJS(multi.effiRE)
multi.effiRE.p = multi.effiRE.p[,c(1:10,110)]
multi.effbssniRE.p = getprobsCJS(multi.effbssniRE)
multi.effbssniRE.p = multi.effbssniRE.p[,c(1:10,110)]
## save
save(multi.dot, multi.dot.p, multi.iRE,  multi.iRE.p, multi.effiRE, multi.effiRE.p, 
     multi.effbssniRE, multi.effbssniRE.p, multimodel.mm4, file="multimodel.npc.cjs.rdata")

# fits when include "precaptures"
## data setup
CH = data.matrix(read.csv("ch.yy.multi.pc.20190315.csv",colClasses="numeric"))
setup = processdata(CH,data.type="sometimes",covs=eff,known=known$kn)
## fit models with all parameters kept
multi.dot.pc = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1,mod.delta=~1,nchains=3,burnin=25000,iter=500000, thin=250,parms="all")
multi.iRE.pc = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+h,mod.delta=~1,nchains=3,burnin=25000,iter=500000, thin=250,parms="all")
multi.effiRE.pc = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff02+h,mod.delta=~1,nchains=3,burnin=25000,iter=500000, thin=250,parms="all")
multi.effbssniRE.pc = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff.bssn+h,mod.delta=~1,nchains=3,burnin=25000,iter=500000, thin=250,parms="all")
# also explored age effect to account for transience in annual apparent survival, but this resulted in poor differentiation of phi from 1:
multi.effbssniRE.age.pc = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff.bssn+h,mod.phi=~1+age,
                                       mod.delta=~1,nchains=3,burnin=25000,iter=500000, thin=250, parms=c("all"), 
                                       parameters=list(Phi=list(age.bins=c(0,1,11))),right=FALSE)
## compare: pdot, pire, peffire, pbssnire
multimodel.mm4.pc = multimodelCJS(list(multi.dot.pc, multi.iRE.pc, multi.effiRE.pc, multi.effbssniRE.pc),
                                  monparms=c("phi","p"))
multimodel.mm4.pc$rjmcmc = multimodel.mm4.pc$rjmcmc[,c(1:10,110:111)]
gelman.diag(multimodel.mm4.pc$rjmcmc, multivariate=F)
heidel.diag(multimodel.mm4.pc$rjmcmc)
effectiveSize(multimodel.mm4.pc$rjmcmc)["M"]
## get p and phi
multi.iRE.pc.p = getprobsCJS(multi.iRE.pc)
multi.iRE.pc.p = multi.iRE.pc.p[,c(1,110)]
multi.effiRE.pc.p = getprobsCJS(multi.effiRE.pc)
multi.effiRE.pc.p = multi.effiRE.pc.p[,c(1:10,110)]
multi.effbssniRE.pc.p = getprobsCJS(multi.effbssniRE.pc)
multi.effbssniRE.pc.p = multi.effbssniRE.pc.p[,c(1:10,110)]
multi.effbssniRE.age.pc.p = getprobsCJS(multi.effbssniRE.age.pc)
multi.effbssniRE.age.pc.p = multi.effbssniRE.age.pc.p[,c(1:10,108:109)]
## save
save(multi.dot.pc, multi.iRE.pc,  multi.iRE.pc.p, multi.effiRE.pc, multi.effiRE.pc.p, 
     multi.effbssniRE.pc, multi.effbssniRE.pc.p, multi.effbssniRE.age.pc, multi.effbssniRE.age.pc.p, 
     multimodel.mm4.pc, file="multimodel.pc.cjs.rdata")

# fit final model with effbssniRE and longer chains 
multi.effbssniRE.pc.long = multimarkCJS(mms=setup,data.type="sometimes",mod.p=~1+eff.bssn+h,mod.delta=~1,
                                        nchains=3,burnin=25000, iter=1000000, thin=50, 
                                        parms=c("pbeta","phibeta","sigma2_zp","delta","alpha","psi"))
plot(density(as.matrix(multi.effbssniRE.pc.long$mcmc[,2])))
gelman.diag(multi.effbssniRE.pc.long$mcmc)
heidel.diag(multi.effbssniRE.pc.long$mcmc)
rmeanplot(multi.effbssniRE.pc.long$mcmc)
multi.effbssniRE.pc.long.p = getprobsCJS(multi.effbssniRE.pc.long)
multi.effbssniRE.pc.long.p = multi.effbssniRE.pc.long.p[,c(1:10,110)]
summary(multi.effbssniRE.pc.long.p)
save(multi.effbssniRE.pc.long, multi.effbssniRE.pc.long.p, file="multimark.cjs.final.rdata")


# model diagnostics: identifiability
library(MCMCvis)
# rename params so can call:
dimnames(multi.effiRE$mcmc[[1]])[[2]] <- c("pbeta","eff02","phibeta","sigma2zp","alpha","psi","delta")
dimnames(multi.effiRE$mcmc[[2]])[[2]] <- c("pbeta","eff02","phibeta","sigma2zp","alpha","psi","delta")
dimnames(multi.effiRE$mcmc[[3]])[[2]] <- c("pbeta","eff02","phibeta","sigma2zp","alpha","psi","delta")
# not clear how to set up for eff02 or delta priors so omit, is sigma2_zp right?
MCMCtrace(multi.effiRE$mcmc, params=c("pbeta","phibeta","sigma2zp","alpha","psi"), 
          priors=matrix( c( rnorm(11700,0,1), rnorm(11700,0,1), 1/rgamma(11700,1,scale=0.01), 
                            rbeta(11700,1,1), rbeta(11700,1,1)), ncol=5))
