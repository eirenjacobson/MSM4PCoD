################ define functions ###############

# Function to simulate capture-recapture data for the Pradel model - using fixed superpopulation and constant lambda. 
# Warning: Returns CH including unobs individuals.
sim.CH.lam.constant <- function(PHI, P, B, D, Nsuper){
  n.occasions <- dim(P)[2]
  CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = Nsuper)
  # Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (t in 1:n.occasions){
    ent.occ <- c(ent.occ, rep(t, B[t]))
  }
  # Simulate survival
  for (i in 1:Nsuper){
    CH.sur[i, ent.occ[i]:n.occasions] <- 1 # Write 1s starting when ind. enters the pop.
  } #i
  for (t in 1:(n.occasions-1)) {
    id = which(as.logical(CH.sur[,t]))
    did = sample(id,D[t]) # let predetermined number of currently living animals die
    CH.sur[did,(t+1):n.occasions] = 0
  } #t
  # Simulate capture
  for (i in 1:Nsuper){
    CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
  } #i
  # Full capture-recapture matrix
  CH <- CH.sur * CH.p
  ## Remove individuals never captured
  #cap.sum <- rowSums(CH)
  #never <- which(cap.sum == 0)
  #CH <- CH[-never,]
  Nt <- colSums(CH.sur) # Actual population size
  return(list(CH=CH, Nt=Nt))
}

#########################
# Function to set population trajectory based on life history inputs
LHpar <- function(N1,lam,nocc,phi) {
  Nmeans = N1*lam^(0:(nocc-1))
  Dmeans = round(Nmeans[1:(nocc-1)])*(1-phi)
  Bmeans = round(Nmeans[2:nocc])-(round(Nmeans[1:(nocc-1)])*phi)
  Nsuper = N1 + sum(round(Bmeans))
  PHI <- matrix(phi, ncol = nocc-1, nrow = Nsuper)
  return(list(Nx=Nmeans,Dx=Dmeans,Bx=Bmeans,Nsuper=Nsuper,PHI=PHI))
}

#############################
# Function to create capture probability matrix
obspar <- function(mean.p, tfe=rep(0,nocc), tre.s2=0, ire.s2=0, nocc, Nsuper){
  logit.iRE <- rnorm(Nsuper,0,ire.s2^0.5)
  logit.p.t <- qlogis(mean.p) + tfe + rnorm(nocc,0,tre.s2^0.5)
  P <- plogis(matrix(logit.p.t, Nsuper, nocc, byrow=T)+logit.iRE)
  return(P)
}

########################
# Function to run set of simulations with constant lambda
sims.lam.constant <- function(skey, nsim,N1,lam,phi,nocc,mean.p,ptFE.s2=0,ptRE.s2=0,piRE.s2=0) {
  sp = LHpar(N1,lam,nocc,phi)
  s = list(ins=list(skey=skey, nsim=nsim, N1=N1, lam=lam, phi=phi,nocc=nocc,mean.p=mean.p,ptFE.s2=ptFE.s2,ptRE.s2=ptRE.s2,piRE.s2=piRE.s2))
  s$ch = array(NA, dim=c(sp$Nsuper,nocc,nsim))
  s$ptfe = matrix(NA,nocc,nsim)
  for (i in 1:nsim) {
    ptFE = rnorm(nocc,0,ptFE.s2^0.5)
    P = obspar(mean.p=mean.p, tfe=ptFE, tre.s2=ptRE.s2, ire.s2=piRE.s2, nocc=nocc, Nsuper=sp$Nsuper)
    temp = sim.CH.lam.constant(PHI=sp$PHI,P=P,B=c(N1,round(sp$Bx)),D=round(sp$Dx),sp$Nsuper)
    s$ch[,,i] = temp$CH
    s$ptfe[,i] = ptFE
  }
  return(s)
}



########################## Run Simulations #########################################

# shared inputs
nsim = 100
N1 = 100
p = 0.085

# SCENARIOS

# Scenario 0: no decline, current time length and effort, tFE+tRE=tvar
s0 = sims.lam.constant(skey="00: lam=1, current effort, tFE+tRE", nsim=nsim, N1=N1, lam = 1,phi = 0.98,nocc = 11,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 1: slow decline, current time length and effort, tFE+tRE=tvar
s1 = sims.lam.constant(skey="01: lam=0.966, current effort, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 11,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 2: fast decline, current time length and effort, tFE+tRE=tvar
s2 = sims.lam.constant(skey="02: lam=0.933, current effort, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 11,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
##### varying time and effort
# Scenario 3: slow decline, double time length, current effort, tFE+tRE=tvar
s3 = sims.lam.constant(skey="03: lam=0.966, 21 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 21,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 4: fast decline, double time length, current effort, tFE+tRE=tvar
s4 = sims.lam.constant(skey="04: lam=0.933, 21 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 21,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 5: slow decline, current time length, double effort, tFE+tRE=tvar
p.x2 = 1.4*p
s5 = sims.lam.constant(skey="05: lam=0.966, 2xeff, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 11,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 6: fast decline, current time length, double effort, tFE+tRE=tvar
s6 = sims.lam.constant(skey="06: lam=0.933, 2xeff, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 11,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)

##### extra time-effort combos to get curves for eff and 2xeff
# Scenario 7: slow decline, third time length and current effort, tFE+tRE=tvar
s7 = sims.lam.constant(skey="07: lam=0.966, 16 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 16,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 8: fast decline, third time length and current effort, tFE+tRE=tvar
s8 = sims.lam.constant(skey="08: lam=0.933, 16 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 16,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 9: slow decline, double time, double effort, tFE+tRE=tvar
s9 = sims.lam.constant(skey="09: lam=0.966, 2xeff, 21 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 21,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 10: fast decline, double time, double effort, tFE+tRE=tvar
s10 = sims.lam.constant(skey="10: lam=0.933, 2xeff, 21 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 21,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 11: slow decline, third time, double effort, tFE+tRE=tvar
s11 = sims.lam.constant(skey="11: lam=0.966, 2xeff, 16 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 16,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)
# Scenario 12: fast decline, third time, double effort, tFE+tRE=tvar
s12 = sims.lam.constant(skey="12: lam=0.933, 2xeff, 16 occ, tFE+tRE", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 16,mean.p = p.x2,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)

##### scenarios with only ptRE, no ptfe 
# importance of covariate
# Scenario 13: slow decline, current time length and effort, no tFE (tRE=tvar)
s13 = sims.lam.constant(skey="13: lam=0.966, current effort, tRE=tvar", nsim=nsim, N1=N1, lam = 0.966,phi = 0.95,nocc = 11,mean.p = p,ptRE.s2=0.35,piRE.s2=0.25)
# Scenario 14: fast decline, current time length and effort, no tFE
s14 = sims.lam.constant(skey="14: lam=0.933, current effort, tRE=tvar", nsim=nsim, N1=N1, lam = 0.933,phi = 0.91,nocc = 11,mean.p = p,ptRE.s2=0.35,piRE.s2=0.25)

# Scenario 15: Halving time of 5 years, current time length and effort, tFE+tRE=tvar
s15 = sims.lam.constant(skey="015: lam=0.87, current effort, tFE+tRE", nsim=nsim, N1=N1, lam = 0.87,phi = 0.85,nocc = 11,mean.p = p,ptFE.s2=0.20,ptRE.s2=0.15,piRE.s2=0.25)

rm(LHpar, sim.CH.lam.constant, sims.lam.constant, 
   obspar, p.x2, N1) #, LHpar.occdesign, sims.lam.constant.occdesign)

save(list=ls(), file="sim.Pradel.scenarios.rdata")



