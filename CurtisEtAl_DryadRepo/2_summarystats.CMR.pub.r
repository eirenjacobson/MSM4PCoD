source("1_ziphius.process.pub.r")

indiv.q30 <- indiv.q %>% filter(Q30.i)
indiv.q30.yy <- indiv.q30 %>% group_by(ID, occ.yy, AgeClass, Sex, FinDmg) %>% 
  summarize(n=n(), mmax = max(c(mR[!is.na(mR)&qR<=2], mL[!is.na(mL)&qL<=2]))) %>% ungroup()


# group size
mean(sight$grpsize)   # 3.0
sd(sight$grpsize)   # 1.8 sd
sd(sight$grpsize)/sqrt(nrow(sight))  # 0.2 se


# association between distinctiveness and capture frequency?
## get highest mark per individual at annual occasion level 
df.idsummary = indiv.q30.yy %>% group_by(ID, Sex, FinDmg) %>% summarize(xm=mean(mmax), ns=n())
with(df.idsummary,plot(xm, jitter(ns)))
df.idsummary %>% mutate(m=round(xm)) %>% group_by(m) %>% summarize(xs=mean(ns), xf=mean(as.numeric(FinDmg)), nx=n())
# m    xs    nx
# 1  1.67     6
# 2  1.29    41
# 3  1.71    35
# 4  1.32    19
# 5  1.58    36


# capture frequencies with precaptures (Table 1)
# create new variables for capture types
indiv.multi.pc %<>%  mutate(T1 = LQ2 & !RQ2, T2 = RQ2 & !LQ2, T4 = RQ2 & LQ2)
## collapse to annual occasions and create CH
indiv.multi.pc.yy <- indiv.multi.pc %>% group_by(ID, occ.yy, linkrl) %>% 
  summarize(T1=sum(T1)>0, T2=sum(T2)>0, T4=sum(T4)>0) %>% 
  mutate(captype = if_else(T4, T4*4, T1*1 + T2*2))
known <- indiv.multi.pc.yy %>% group_by(ID) %>% summarize(kn=any(captype==4 | linkrl))
indiv.multi.pc.yy %<>%  arrange(occ.yy,ID)
CH.multi.pc <- spread(indiv.multi.pc.yy[c("ID","occ.yy","captype")], key=occ.yy, value=captype, fill=0, drop=FALSE) %>% 
  ungroup() %>% select(-ID)
## summaries
a <- data.frame(t(table(rowSums(CH.multi.pc[known$kn,]>0))))
#1  2  3   4  5  7
#54 16 11  4  1  1
with(a, data.frame(sum(Freq), sum(as.integer(levels(Var2))*Freq)/sum(Freq)))  #
sides <- indiv.multi.pc.yy %>% group_by(ID) %>% 
  summarize(r=all(captype==2 & !linkrl), l=all(captype==1 & !linkrl))
a <- data.frame(t(table(rowSums(CH.multi.pc[sides$r,]>0))))
#1  2 
#18  3 
with(a, data.frame(sum(Freq), sum(as.integer(Var2)*Freq)/sum(Freq)))  #
a <- data.frame(t(table(rowSums(CH.multi.pc[sides$l,]>0))))
# 1  2 
#17  1
with(a, data.frame(sum(Freq), sum(as.integer(Var2)*Freq)/sum(Freq)))  #

# capture history statistics
df.sex <- indiv.q30 %>% mutate(occ.yy=year(occ.yy)) %>% group_by(ID, Sex) %>% 
  summarise(nocc = n_distinct(occ.yy), 
            nsight = n_distinct(SightID),
            obsdur = max(occ.yy)-min(occ.yy) + 1, 
            meanint = ifelse(n_distinct(occ.yy) > 1, mean(diff(unique(occ.yy))), NA)) %>% 
  ungroup()
## M vs F
df.sex %>% group_by(Sex) %>% summarize(nocc.x = mean(nocc), nocc.se = sd(nocc)/sqrt(n()), 
                                       nsight.x = mean(nsight), nsight.se = sd(nsight)/sqrt(n()),
                                       obsdur.x = mean(obsdur), obsdur.se = sd(obsdur)/sqrt(n()),
                                       meanint.x = mean(meanint, na.rm = T), meanint.se = sd(meanint, na.rm = T)/sqrt(n()))
## overall statistics
table(df.sex$obsdur[df.sex$nocc>1])
#  2  3  4  5  6  7  8  9 10 11 
#  6  6  3  7  2  3  4  3  4  2 
mean(df.sex$meanint[df.sex$nocc>1])  # 3.1
rm(df.sex)

# mean annual sex proportions 
round(colMeans(t(apply(table(indiv.q30.yy$occ.yy,indiv.q30.yy$Sex),1,function(x) x/sum(x)))),3)
#Female    Male Unknown 
#0.453   0.385   0.162 

# mean annual Age Class proportions
round(colMeans(t(apply(table(indiv.q30.yy$occ.yy,indiv.q30.yy$AgeClass),1,function(x) x/sum(x)))),3)
#Adult              Calf          Juvenile          Subadult           Unknown   Unknown, immature 
#0.615            0.047             0.085             0.215              0.035          0.004 

# IDs of Age Classes by year for included records
temp <- indiv.sel %>% filter(Q30.i) %>% 
  group_by(ID, occ.yy, ssn, AgeClass, Sex, FinDmg) %>% 
  summarize(n=n(), mmax = max(c(mR[!is.na(mR)&qR<=2], mL[!is.na(mL)&qL<=2]))) %>% ungroup()
round(colMeans(t(apply(table(temp$occ.yy,temp$AgeClass),1,function(x) x/sum(x)))),3)
#Adult          Juvenile          Subadult           Unknown Unknown, immature 
#0.659             0.079             0.221             0.036             0.005 

# calf accompaniment
# rejected:
#(1) num calves/num yrs sighted for identified females (Weinrich and Corbelli 2009) 
#    - don't have these data on hand (is this association even reliably apparent to data collectors?)
#(2) Barlow 1990 (?): birth interval approach: use all females with previous births to estimate birth interval
#    - same problem as (1)
#(3) reproductive rate: # calves with known females/# known females 
#    - this may not be plausible since may not know adult female with high confidence until has had calf.
#    - also see above
#(4) Steiger and Calambokidis 2000:
#   (b) repro rate based on photo-identification: % mothers (based on calf) id'd/ all whales id'd/ yr
#     - see above issues
#   (c) max calving rate: calves per mature female per year, using CHs of females from first year seen with calf (Clapham and Mayo 1990)
#     - see above issues
#   (d) alternate calving rate: max calving rate is biased high (Barlow and Clapham 1997) so start with year after first year seen with calf - see Barlow 1990 above

#(1) rate of calf production: divide num calves by total indiv in data set (Ratnaswamy and Winn 1993)
#  -> see annual level summary above, which accounts for individual transitions in age class

#(2) Steiger and Calambokidis 2000: 
#   (a) repro rate based on sightings:  proportion of all whales approached that were calves (use photos or sightings)   
repro.rate <- indiv.q %>% group_by(AgeClass) %>% summarise(n=n()) %>% ungroup() %>% mutate(p=round(n/sum(n),3))
binom.confint(repro.rate$n[repro.rate$AgeClass=="Calf"],sum(repro.rate$n),conf.level=0.95,methods="wilson")
#AgeClass              n     p
# Adult               170 0.637
# Calf                 16 0.06 
# Juvenile             21 0.079
# Subadult             44 0.165
# Unknown              15 0.056
# Unknown, immature     1 0.004

