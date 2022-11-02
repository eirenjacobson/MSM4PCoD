library(rjags)
library(mcmcplots)

### read data
CH.R = data.matrix(read.csv("ch.yy.Prad.R.csv",colClasses="numeric"))
occ.cov = read.csv("occasions.cov.yy.20190315.csv")

### derive data for the model
nocc <- ncol(CH.R)
# e is index of the earliest observation
get.first <- function (x) min( which (x !=0))
e.R <- apply (CH.R ,1, get.first)
# l is index of the last observation
get.last <- function (x) max( which (x !=0))
l.R <- apply (CH.R ,1 ,get.last)
# u = number of animals observed for the first time at i
u.R <- hist(e.R,breaks=0:nocc+0.5, plot=F)$counts
# n = number of animals observed at i
n.R <- colSums(CH.R)
# v = number of animals observed for the last time at i
v.R <- hist(l.R,breaks=0:nocc+0.5, plot=F)$counts
# d = number of animals removed from the population at time i
d <- rep(0, nocc)

## ptRE
eff <- rep(0,11)
# Bundle data
bugs.data.R <- list (u=u.R, n=n.R, v=v.R, d=d, s=nocc, eff = as.vector(eff))
## Initial values
inits <- function (){ list (
  beta.eff.p = runif (1, -2, 2),
  alpha.rho = runif (1, -0.5, 0.5) ,
  mu = runif (nocc, 0.3 , 1)
)}
# Define parameters to be monitored
parameters <- c("mean.phi", "alpha.phi",
                "gamma[2]" ,
                "alpha.rho", "mean.rho", 
                "p", "mean.p","alpha.p", #"beta.eff.p", 
                "sigma.p")
# MCMC settings
ni <-32500; nt <- 1; nc <- 3; nad <- 500; nb <- 2000
# Call JAGS
# right side
m <- jags.model(file="PradelLambda.pefftRE.alt.bug", data=bugs.data.R, inits=inits, n.chains=nc, n.adapt=nad)
update(m, nb)
draws.PL.ptRE.R <- coda.samples(m, ni, variable.names=parameters, thin=nt)
gelman.diag(draws.PL.ptRE.R,multivariate=F)
heidel.diag(draws.PL.ptRE.R)
plot(draws.PL.ptRE.R[,c("mean.phi","mean.p","mean.rho")])
rmeanplot(draws.PL.ptRE.R)
summary(draws.PL.ptRE.R)

## pefftRE
eff <- occ.cov$eff.bssn
# Bundle data
bugs.data.R <- list (u=u.R, n=n.R, v=v.R, d=d, s=nocc, eff = as.vector(eff))
## Initial values
inits <- function (){ list ( #alpha.phi = runif (1, -2, 2),
  beta.eff.p = runif (1, -2, 2),
  alpha.rho = runif (1, -0.5, 0.5) ,
  #alpha.p = runif (1, -0.5, 0.5) ,
  mu = runif (nocc, 0.3 , 1)
)}
# Define parameters to be monitored
parameters <- c("mean.phi", "alpha.phi",
                "gamma[2]" ,
                "alpha.rho", "mean.rho", 
                "p", "mean.p","alpha.p","beta.eff.p", 
                "sigma.p")
# MCMC settings
ni <-32500; nt <- 1; nc <- 3; nad <- 500; nb <- 2000
# Call JAGS
m <- jags.model(file="PradelLambda.pefftRE.alt.bug", data=bugs.data.R, inits=inits, n.chains=nc, n.adapt=nad)
update(m, nb)
draws.PL.pefftRE.R <- coda.samples(m, ni, variable.names=parameters, thin=nt)
gelman.diag(draws.PL.pefftRE.R,multivariate=F)
heidel.diag(draws.PL.pefftRE.R)
plot(draws.PL.pefftRE.R[,c("mean.phi","mean.p","mean.rho")])
rmeanplot(draws.PL.pefftRE.R)
summary(draws.PL.pefftRE.R)


## pefftRE phifix: multiple imputation approach
eff <- occ.cov$eff.bssn
## Initial values
inits <- function (){ list ( beta.eff.p = runif (1, -2, 2),
                             alpha.rho = runif (1, -0.5, 0.5) ,
                             mu = runif (nocc, 0.3 , 1)
)}
# Define parameters to be monitored
parameters <- c("gamma[2]" , "mean.phi",
                "alpha.rho", "mean.rho", 
                "p", "alpha.p", "mean.p","beta.eff.p", 
                "sigma.p")
# MCMC settings
ni <- 20000; nt <- 100; nc <- 3; nad <- 500;  nb <- 10000
ndi = 250
attach("multimark.cjs.final.rdata")
phifix = sample(as.matrix(multi.effbssniRE.pc.long.p[,"phi[10,10]"]), size=ndi)
detach(2)
matdraws.PL.pefftRE.phifix <- as.data.frame(matrix(NA, ndi*600, 19))   
for (di in 1:ndi) {
  # Bundle data
  bugs.data <- list (u=u.R, n=n.R, v=v.R, d=d, s=nocc, eff = as.vector(eff), phifix=phifix[di])
  # Call JAGS
  m <- jags.model(file="PradelLambda.pefftRE.fixedphi.bug", data=bugs.data, inits=inits, n.chains=nc, n.adapt=nad)
  update(m, nb)
  temp <- coda.samples(m, ni, variable.names=parameters, thin=nt)
  matdraws.PL.pefftRE.phifix[(di*600-599):(di*600),] <- as.matrix(temp)
}
names(matdraws.PL.pefftRE.phifix) <- dimnames(temp[[1]])[[2]]
draws.PL.pefftRE.phifix.R <- as.mcmc.list(mcmc(matdraws.PL.pefftRE.phifix))
# iteratively plot rho (lambda) distributions to check convergence
plot(density(as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"])[1:15000]),ylim=c(0,14))
for (i in 2:9) {
  lines(density(as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"])[1:(i*15000)]),col=3)
}
lines(density(as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"])[1:(10*15000)]),col=3)
rmeanplot(draws.PL.pefftRE.phifix.R)
summary(as.mcmc.list(mcmc(draws.PL.pefftRE.phifix.R[[1]][1:100000,])))
summary(draws.PL.pefftRE.phifix.R)


save(draws.PL.ptRE.R, draws.PL.ptRE.L, draws.PL.pefftRE.R, draws.PL.pefftRE.L, 
     draws.PL.pefftRE.phifix.R, file="PradelLambda.20190315.rdata")
