# Bayesian Pradel-lambda analyses for lambda simulations

attach("sim.Pradel.scenarios.20190315.rdata")

library(rjags)
library(dplyr)
library(readr)

# define data prep functions
get.first <- function (x) min( which (x !=0))
get.last <- function (x) max( which (x !=0))
# MCMC settings
nit <- 10000; nc <- 3; nad <- 500; nb <- 1000;
# Define parameters to be monitored
parameters <- c("mean.phi", "mean.rho", "mean.p") #,"sigma.p")
## Initial values
inits <- function (){ list ( 
  alpha.phi = runif (1, 0, 5),
  alpha.p = runif (1, -5, 0),
  beta.eff.p = runif (1, -2, 2),
  alpha.rho = runif (1, -0.5, 0.5),
  mu = runif (ncol(ch) , 0.3 , 1)
)}

scenarios = 0:15
models.summary = data.frame(NULL)

for (s in scenarios) {
  for (i in 1:nsim) {
    
    ## read data
    chall <- get(paste("s",s,sep=""))$ch[,,i]
    ch <- chall[-which(rowSums(chall)==0),]
    if (get(paste("s",s,sep=""))$ins$ptFE.s2 > 0) {
      eff <- get(paste("s",s,sep=""))$ptfe[,i]
    } else { 
      eff <- rep(0, ncol(ch)) 
    }
    
    ### derive data for the model
    # e is index of the earliest observation
    e <- apply (ch ,1, get.first)
    # l is index of the last observation  
    l <- apply (ch ,1 ,get.last)
    # u = number of animals observed for the first time at i
    u <- hist(e,breaks=0:ncol(ch)+0.5,plot=F)$counts
    # n = number of animals observed at i
    n <- colSums(ch)
    # v = number of animals observed for the last time at i
    v <- hist(l,breaks=0:ncol(ch)+0.5,plot=F)$counts
    # d = number of animals removed from the population at time i
    d <- rep(0, dim(ch)[2])
    # Bundle data
    bugs.data <- list (u=u, n=n, v=v, d=d, s=ncol(ch), eff = as.vector(eff))
    
    # Call JAGS
    m <- jags.model(file="PradelLambda.pefftRE.alt.bug", data=bugs.data, inits=inits, n.chains=nc, n.adapt=nad)
    update(m, nb)
    draws <- coda.samples(m, nit, variable.names=parameters)
    gd <- gelman.diag(draws, multivariate=F)
    matdraws <- data.frame(as.matrix(draws))
    
    for (par in c("phi","p","rho")) {
      x <- matdraws[,paste("mean",par,sep=".")]
      den <- density(x, kernel=c("gaussian"),from=min(x),to=max(x))
      denranks = order(den$y,decreasing=T)
      cum.perc.ranked = cumsum(den$y[denranks])/sum(den$y)
      assign(paste(par,"mode",sep="."), den$x[den$y==max(den$y)])
      assign(paste(par,"hpdi.90",sep="."), range(den$x[denranks][!(cum.perc.ranked>0.9)]))
      assign(paste(par,"hpdi.95",sep="."), range(den$x[denranks][!(cum.perc.ranked>0.95)]))
    }
    
    # pull together results: 
    tempsumm <- with(matdraws, 
      data.frame(scenario=s, sim=i, psrf.max=max(gd$psrf),
                phi.x=mean(mean.phi), phi.se=sd(mean.phi), phi.med=median(mean.phi), 
                phi.lcl=quantile(mean.phi,0.025), phi.ucl=quantile(mean.phi,0.975), 
                phi.10=quantile(mean.phi,0.1), phi.90=quantile(mean.phi,0.9), 
                phi.05=quantile(mean.phi,0.05), phi.95=quantile(mean.phi,0.95),
                phi.mode=phi.mode, phi.den.05=phi.hpdi.90[1], phi.den.95=phi.hpdi.90[2],
                phi.den.025=phi.hpdi.95[1], phi.den.975=phi.hpdi.95[2], 
                p.x=mean(mean.p), p.se=sd(mean.p), p.med=median(mean.p), 
                p.lcl=quantile(mean.p,0.025), p.ucl=quantile(mean.p,0.975), 
                p.10=quantile(mean.p,0.1), p.90=quantile(mean.p,0.9), 
                p.05=quantile(mean.p,0.05), p.95=quantile(mean.p,0.95), 
                p.mode=p.mode, p.den.05=p.hpdi.90[1], p.den.95=p.hpdi.90[2], 
                p.den.025=p.hpdi.95[1], p.den.975=p.hpdi.95[2], 
                lam.x=mean(mean.rho), lam.se=sd(mean.rho), lam.med=median(mean.rho), 
                lam.lcl=quantile(mean.rho,0.025), lam.ucl=quantile(mean.rho,0.975), 
                lam.10=quantile(mean.rho,0.1), lam.90=quantile(mean.rho,0.9), 
                lam.05=quantile(mean.rho,0.05), lam.95=quantile(mean.rho,0.95), 
                lam.mode=rho.mode, lam.den.05=rho.hpdi.90[1], lam.den.95=rho.hpdi.90[2],
                lam.den.025=rho.hpdi.95[1], lam.den.975=rho.hpdi.95[2], 
                lam.pdec=sum(mean.rho<1)/nrow(matdraws)))
                
    models.summary = plyr::rbind.fill(models.summary, tempsumm)
    rm(chall,ch,eff,e,l,u,n,v,d,bugs.data,m,draws,gd,matdraws,par,x,den,denranks,cum.perc.ranked,tempsumm,
       phi.mode, phi.hpdi.90, phi.hpdi.95, p.mode, p.hpdi.90, p.hpdi.95, 
       rho.mode, rho.hpdi.90, rho.hpdi.95)
  }
  ni = (s*nsim+1):((s+1)*nsim)
  
  # add true values to data frame
  ins = get(paste("s",s,sep=""))$ins
  models.summary$skey[ni] = ins$skey
  models.summary$lam.true[ni] = ins$lam
  models.summary$phi.true[ni] = ins$phi
  models.summary$p.true[ni] = ins$mean.p
  models.summary$nocc[ni] = ins$nocc
  models.summary$ny[ni] = ifelse(is.null(ins$ny), ins$nocc, ins$ny)
  models.summary$ptfe.s2[ni] = ins$ptFE.s2
  models.summary$ptre.s2[ni] = ins$ptRE.s2
  models.summary$pire.s2[ni] = ins$piRE.s2
  rm(ins)
  
  if(s==0) {
    write_csv(models.summary[ni,], path="sims.lam.Bayes.Pradel.tRE.csv")
  } else{
    write_csv(models.summary[ni,], path="sims.lam.Bayes.Pradel.tRE.csv",append=T)
  }
}
rm(i, s, ni, nit, nt, nc, nad, parameters, get.first, get.last, scenarios)
 
save(list=ls(), file="sims.lam.Bayes.Pradel.tRE.rdata")