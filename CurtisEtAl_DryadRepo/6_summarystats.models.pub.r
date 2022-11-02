# assemble data frame with mode, HPDI at 75% and 90%, median, 75 and 90% percentiles, mean, sd

####################### DEFINE FUNCTIONS #############################

# define  function to extract stats
summstats <- function(x) {
  # calculate mode and hpdi
  den <- density(x, kernel=c("gaussian"),from=min(x),to=max(x))
  mode = den$x[den$y==max(den$y)]
  denranks = order(den$y,decreasing=T)
  cum.perc.ranked = cumsum(den$y[denranks])/sum(den$y)
  hpdi.50 = range(den$x[denranks][!(cum.perc.ranked>0.50)])
  hpdi.90 = range(den$x[denranks][!(cum.perc.ranked>0.9)])
  # calculate median and percentiles
  percs = as.numeric(quantile(x,probs=c(0.05,.25,0.5,.75,0.95)))
  return(round(data.frame(hpdi.90.l=hpdi.90[1],hpdi.50.l=hpdi.50[1],mode=mode,hpdi.50.u=hpdi.50[2],hpdi.90.u=hpdi.90[2],
              perc.05=percs[1],perc.25=percs[2],perc.50=percs[3],perc.75=percs[4],perc.95=percs[5],
              mean=mean(x), sd=sd(x)),3))
}

# process mcmc.list object 
process.mcmc.list.stats <- function(x) {
  all = as.matrix(x)
  nv = dim(all)[2]
  summarystats = data.frame(NULL)
  for (n in 1:nv) {
    summarystats = plyr::rbind.fill(summarystats, summstats(all[,n]))
    summarystats$par[n] = dimnames(all)[[2]][n]
  }
  summarystats = summarystats[c("par",names(summarystats)[1:(dim(summarystats)[2]-1)])]
  return(summarystats)
}

  
#####################################################################

# load packages
library(rjags)
library(dplyr)
library(tidyr)

# summary stats for multimark model choice
attach("multimodel.pc.cjs.rdata")
out = data.frame(model=names(multimodel.mm4.pc$pos.prob$overall), 
                 prob=round(as.vector(multimodel.mm4.pc$pos.prob$overall),3))
neff = round(effectiveSize(multimodel.mm4.pc$rjmcmc)["M"])
gd.mv = round(gelman.diag(multimodel.mm4.pc$rjmcmc)$mpsrf,3)
out = out %>% spread(key=model, value=prob)
out = cbind(out, gd.mv, neff)
write.csv(out, file="summarystats.multimark.cjs.multimodel.20190315.csv",row.names=F)
detach(2)

# summary stats for multimark CJS model
attach("multimark.cjs.final.rdata")
out = process.mcmc.list.stats(multi.effbssniRE.pc.long$mcmc)
out.p = process.mcmc.list.stats(multi.effbssniRE.pc.long.p)
out <- rbind(out, out.p)
neff = round(c( effectiveSize(multi.effbssniRE.pc.long$mcmc),
          effectiveSize(multi.effbssniRE.pc.long.p)))
out <- cbind(out, neff)
write.csv(out, file="summarystats.multimark.cjs.20190315.csv",row.names=F)
detach(2)

# summary stats for multimark closed-population model
attach("multimark.closed.20190315.rdata")
attach("cf.rdata")   # correction factor from ziphius.selectdata.r
out = process.mcmc.list.stats(multi.effiRE.sub$mcmc)
out.p = process.mcmc.list.stats(multi.effiRE.sub.p)
out <- rbind(out, out.p)
out[3,2:ncol(out)] <- cf.nc*out[3,2:ncol(out)]   
neff = round(c( effectiveSize(multi.effiRE.sub$mcmc),
          effectiveSize(multi.effiRE.sub.p)))
out <- cbind(out, neff)
write.csv(out, file="summarystats.multimark.closed.cf.nc.csv",row.names=F)
out = process.mcmc.list.stats(multi.effiRE.sub$mcmc)
out.p = process.mcmc.list.stats(multi.effiRE.sub.p)
out <- rbind(out, out.p)
out[3,2:ncol(out)] <- cf.c*out[3,2:ncol(out)]   
neff = round(c( effectiveSize(multi.effiRE.sub$mcmc),
          effectiveSize(multi.effiRE.sub.p)))
out <- cbind(out, neff)
write.csv(out, file="summarystats.multimark.closed.cf.c.csv",row.names=F)
detach(2)
detach(2)

# summary stats for Pradel-lambda model
attach("PradelLambda.20190315.rdata")
# pefftRE phifix
out = process.mcmc.list.stats(draws.PL.pefftRE.phifix.R)
#pr.neg.phifix <- sum(as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"])<1)/length(as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"]))   # 0.624
neff = round(effectiveSize(draws.PL.pefftRE.phifix.R))
out <- cbind(out, neff)
write.csv(out, file="summarystats.R.PradelLambda.pefftRE.phifix.csv",row.names=F)
# pefftRE
out = process.mcmc.list.stats(draws.PL.pefftRE.R)
#pr.neg.pefftRE <- sum(as.matrix(draws.PL.pefftRE.R[,"mean.rho"])<1)/length(as.matrix(draws.PL.pefftRE.R[,"mean.rho"]))   # 0.719
neff = round(effectiveSize(draws.PL.pefftRE.R))
out <- cbind(out, neff)
write.csv(out, file="summarystats.R.PradelLambda.pefftRE.csv",row.names=F)
# ptRE
out = process.mcmc.list.stats(draws.PL.ptRE.R)
neff = round(effectiveSize(draws.PL.ptRE.R))
out <- cbind(out, neff)
#pr.neg.ptRE <- sum(as.matrix(draws.PL.ptRE.R[,"mean.rho"])<1)/length(as.matrix(draws.PL.ptRE.R[,"mean.rho"]))   # 0.307
write.csv(out, file="summarystats.R.PradelLambda.ptRE.csv",row.names=F)
detach(2)

