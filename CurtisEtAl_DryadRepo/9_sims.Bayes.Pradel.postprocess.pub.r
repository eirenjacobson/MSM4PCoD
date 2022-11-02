# post-processsing of Bayesian Pradel model sims

library(readr)
library(dplyr)
library(tidyr)
library(magrittr)

### import data (either (a) OR (b)
## (a) combine results files
ct <- cols(
  scenario = col_integer(),
  sim = col_integer(),
  psrf.max = col_double(),
  phi.x = col_double(),
  phi.se = col_double(),
  phi.med = col_double(),
  phi.lcl = col_double(),
  phi.ucl = col_double(),
  phi.10 = col_double(),
  phi.90 = col_double(),
  phi.05 = col_double(),
  phi.95 = col_double(),
  phi.mode = col_double(),
  phi.den.05 = col_double(),
  phi.den.95 = col_double(),
  phi.den.025 = col_double(),
  phi.den.975 = col_double(),
  p.x = col_double(),
  p.se = col_double(),
  p.med = col_double(),
  p.lcl = col_double(),
  p.ucl = col_double(),
  p.10 = col_double(),
  p.90 = col_double(),
  p.05 = col_double(),
  p.95 = col_double(),
  p.mode = col_double(),
  p.den.05 = col_double(),
  p.den.95 = col_double(),
  p.den.025 = col_double(),
  p.den.975 = col_double(),
  lam.x = col_double(),
  lam.se = col_double(),
  lam.med = col_double(),
  lam.lcl = col_double(),
  lam.ucl = col_double(),
  lam.10 = col_double(),
  lam.90 = col_double(),
  lam.05 = col_double(),
  lam.95 = col_double(),
  lam.mode = col_double(),
  lam.den.05 = col_double(),
  lam.den.95 = col_double(),
  lam.den.025 = col_double(),
  lam.den.975 = col_double(),
  lam.pdec = col_double(),
  skey = col_character(),
  lam.true = col_double(),
  phi.true = col_double(),
  p.true = col_double(),
  nocc = col_integer(),
  ny = col_integer(),
  ptfe.s2 = col_double(),
  ptre.s2 = col_double(),
  pire.s2 = col_double()
)

## (b) read in results
models.summary <- read_csv("sims.lam.Bayes.Pradel.pefftRE.csv", col_types=ct)
rm(ct)

## models converged?
models.summary %<>% mutate(converged=psrf.max<1.05)
mconv <- models.summary %>% group_by(skey) %>% summarize(pconverged = sum(converged)/n())

## summaries
sims.stats <- models.summary %>% filter(converged) %>% 
  group_by(skey) %>% 
  summarize(dec.80 = sum(lam.pdec>=0.8)/n(), 
            dec.95 = sum(lam.pdec>=0.95)/n(),
            dec.CI = sum(lam.ucl<1)/n(),
            dec.hpdi.95 = sum(lam.den.975<1)/n(),
            
            inc.80 = sum(lam.pdec<=0.2)/n(),
            inc.CI = sum(lam.lcl>1)/n(),
            
            lambias.x = mean(lam.x-lam.true), 
            lambias.med = mean(lam.med-lam.true), 
            lambias.mode = mean(lam.mode-lam.true), 
            lam.inCI = sum((lam.lcl<=lam.true) & (lam.ucl>=lam.true))/n(),
            lam.in80 = sum((lam.10<=lam.true) & (lam.90>=lam.true))/n(),
            lam.in90 = sum((lam.05<=lam.true) & (lam.95>=lam.true))/n(),
            lam.inhpdi = sum((lam.den.025<=lam.true) & (lam.den.975>=lam.true))/n(),
            lam.inhpdi.90 = sum((lam.den.05<=lam.true) & (lam.den.95>=lam.true))/n(),
            
            phibias.x = mean(phi.x-phi.true), 
            phibias.med = mean(phi.med-phi.true), 
            phibias.mode = mean(phi.mode-phi.true),
            phi.inCI = sum((phi.lcl<=phi.true) & (phi.ucl>=phi.true))/n(),
            phi.in80 = sum((phi.10<=phi.true) & (phi.90>=phi.true))/n(),
            phi.in90 = sum((phi.05<=phi.true) & (phi.95>=phi.true))/n(),
            phi.inhpdi = sum((phi.den.025<=phi.true) & (phi.den.975>=phi.true))/n(),
            phi.inhpdi.90 = sum((phi.den.05<=phi.true) & (phi.den.95>=phi.true))/n()) %>% 
  cbind(mconv$pconverged)
rm(mconv)

# Calculate inference error curves for current effort scenarios
lam1.error1 <- models.summary %>% filter(converged & scenario==0) %>% select(lam.pdec) %>% 
  summarize(a.50=sum(lam.pdec>=0.5)/n(), a.55=sum(lam.pdec>=0.55)/n(), a.60=sum(lam.pdec>=0.6)/n(), a.65=sum(lam.pdec>=0.65)/n(), 
            a.70=sum(lam.pdec>=0.7)/n(), a.75=sum(lam.pdec>=0.75)/n(), a.80=sum(lam.pdec>=0.8)/n(), a.85=sum(lam.pdec>=0.85)/n(), 
            a.90=sum(lam.pdec>=0.9)/n(), a.95=sum(lam.pdec>=0.95)/n()) %>% 
  pivot_longer(names_prefix="a.", cols = a.50:a.95, values_to="error I", names_to = "a") %>% 
  mutate(a=as.numeric(a)/100)
  
lam.dec.error2 <- models.summary %>% filter(converged & scenario %in% c(1:2,15)) %>% group_by(scenario) %>% select(scenario, lam.pdec) %>% 
  summarize(a.50=1-sum(lam.pdec>=0.5)/n(), a.55=1-sum(lam.pdec>=0.55)/n(), a.60=1-sum(lam.pdec>=0.6)/n(), a.65=1-sum(lam.pdec>=0.65)/n(), 
            a.70=1-sum(lam.pdec>=0.7)/n(), a.75=1-sum(lam.pdec>=0.75)/n(), a.80=1-sum(lam.pdec>=0.8)/n(), a.85=1-sum(lam.pdec>=0.85)/n(), 
            a.90=1-sum(lam.pdec>=0.9)/n(), a.95=1-sum(lam.pdec>=0.95)/n()) %>%  
  pivot_longer(names_prefix="a.", cols = a.50:a.95, values_to="error II", names_to = "a") %>% 
  mutate(a=as.numeric(a)/100) %>% 
  pivot_wider(names_from="scenario", values_from="error II") %>% mutate(a=(10:19)/20)
