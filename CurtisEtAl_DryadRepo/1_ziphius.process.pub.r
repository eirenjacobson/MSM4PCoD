################ PROCESS ZIPHIUS MARK-RECAPTURE DATA #################

## load general packages
#library(tidyverse)
library(readr)
library(dplyr)
library(lubridate)
library(magrittr)
library(tidyr)

## import data
library(readxl)
datapath <- "SCORE-Zc-Data 20190315.xlsx"
sight <- read_excel(datapath, sheet="SCORE-Zc-Sightings", col_names=TRUE)
indiv <- read_excel(datapath, sheet="SCORE-Zc-ID", col_names=TRUE)
detach("package:readxl")
effort <- readr::read_csv("tSCORE-Zica-Effort 20190315.txt", 
                          col_types = cols_only(
                            Vessel = col_character(),
                            Date = col_character(),
                            Time = col_character(),
                            Latitude = col_double(),
                            Longitude = col_double(),
                            Depth = col_integer(),
                            `800m` = col_integer(),
                            CatalinaBasin = col_double(),
                            SanNic = col_double(),
                            SantaCruz = col_double(),
                            SOAR = col_double(),
                            Q = col_character(),
                            B = col_integer(),
                            V = col_double(),
                            S = col_integer()
                          ))

# load old effort file to splice in effort for Vessel N2 10/24/2007 after 11:11
effort.20180427 <- readr::read_csv("SCORE-Zica-TrackGISCond 20180427.txt", 
                                   col_types = cols_only(
                                     Vessel = col_character(),
                                     Date = col_character(),
                                     Time = col_character(),
                                     Latitude = col_double(),
                                     Longitude = col_double(),
                                     Depth = col_integer(),
                                     `800m` = col_integer(),
                                     Cat = col_double(),
                                     SN = col_double(),
                                     SC = col_double(),
                                     SOAR = col_double(),
                                     Q = col_character(),
                                     B = col_integer(),
                                     V = col_double(),
                                     S = col_integer()
                                   ))

## change data field names
sight  %<>%  rename(Time='Start Time', lat='Start Dec Lat', lon='Start Dec Long', grpsize='Grp Size')
indiv %<>% rename(qL='Qual L', qR='Qual R', mL='Marks L', mR='Marks R', q3L='Q3 L', q3R='Q3 R')
effort %<>% rename(X800m='800m', Cat="CatalinaBasin", SN="SanNic", SC="SantaCruz")
effort.20180427 %<>% rename(X800m='800m')

# strip time zone from dates read in from worksheet
sight %<>% mutate(Date=date(Date))
indiv %<>% mutate(Date=as.Date(Date, format="%d-%b-%y"))

## create sighting-level link field
sight %<>% mutate(SightID = paste(Date, Vessel, Sighting, sep="_"))
indiv %<>% mutate(SightID = paste(Date, Vessel, Sighting, sep="_"))

## edit sightings data
# create date-time field
sight %<>% rowwise() %>% mutate(DateTime = as.POSIXct(paste(Date, Time, sep=" "))) %>% ungroup()

## edit effort data
# create date-time and date fields
effort  %<>%  rowwise() %>% 
  mutate(
    DateTime = as.POSIXct(paste(
    sub("0:00:00", "", Date),
    unlist(strsplit(Time, " "))[2],
    sep = ""
    ), format = "%m/%d/%Y %H:%M:%S"),
    Date = date(floor_date(DateTime, "day"))
  ) %>% ungroup()
effort.20180427  %<>%  rowwise() %>% 
  mutate(
    DateTime = as.POSIXct(paste(
      sub("0:00:00", "", Date),
      unlist(strsplit(Time, " "))[2],
      sep = ""
    ), format = "%m/%d/%Y %H:%M:%S"),
    Date = date(floor_date(DateTime, "day"))
  ) %>% ungroup()

# splice in data from effort.20180427, "turn off" effort at last record before that in effort
replacement <- tail(effort[effort$Vessel=="N2" & effort$DateTime<as.POSIXct("2007-10-24 11:11"),],1)
ri <- which(duplicated(bind_rows(replacement, effort)))-1
replacement[,c("X800m","SN")]=0
effort[ri,] <- replacement
effort <- effort.20180427 %<>% filter(Vessel=="N2" & DateTime>as.POSIXct("2007-10-24 11:11") & Date<date("2007-10-25")) %>% 
  bind_rows(effort)
rm(effort.20180427)
rm(replacement,ri)

# calculate minutes of effort for each GPS record
effort %<>% arrange(Vessel, DateTime)
effort %<>% group_by(Vessel, Date) %>%
  mutate(
    Minutes = as.numeric(difftime(lead(DateTime), DateTime, units = "mins"))) %>% 
  ungroup()
# limit effort records to inside SanNic Basin and >800m (all Zc sightings are at >=800m)
effort <- effort %>% dplyr::filter(SN>0 & X800m>0)

### omit sightings data not appropriate for Ziphius CMR analysis
## drop sightings data from outside San Nicolas
sight %<>% dplyr::filter(Basin=="SanNic")
indiv %<>% dplyr::filter(SightID %in% sight$SightID)
## drop sighting for 2011-05-05 (Tanner Cyn, adjacent to but not audible from SOAR array)
# No individual IDs from that sighting, so only edit group data; effort data coded as !SN
sight %<>% dplyr::filter(SightID!="2011-05-05_PHY_5")

## add effort type to sightings
sight <- mutate(sight, Q = NA_character_, B = NA_integer_, V = NA_real_, S = NA_integer_, 
                Depth = NA_real_, efflon = NA_real_, efflat = NA_real_)
effint <- with(effort, interval(DateTime,DateTime+dminutes(Minutes)-dseconds(1)))
for (si in 1:nrow(sight)) {
  ei <- which((sight$DateTime[si] %within% effint) & sight$Vessel[si]==effort$Vessel)
  if (length(ei)) sight[si,c("Q","B","V","S","Depth","efflon","efflat")] <- effort[ei,c("Q","B","V","S","Depth","Longitude","Latitude")]
}
rm(effint,si,ei)

## further data selection
# omit data before 2007 captures (before 2007-Oct), after 2018-Jul 
# (eliminates low-sightings occasions in March 2006, April 2007 at beginning of time series that could drive pent fit in POPAN models)
sight <- dplyr::filter(sight, Date>=as.POSIXct("2007-08-01") & Date<as.POSIXct("2018-08-01"))
indiv %<>% dplyr::filter(Date>=as.POSIXct("2007-08-01") & Date<as.POSIXct("2018-08-01"))
effort  <- dplyr::filter(effort, Date>=as.POSIXct("2007-08-01") & Date<as.POSIXct("2018-08-01"))

## set seasons, occasions, add to data frames
ssn <- 7:10
occasions.yy <- date(unique(floor_date(effort$Date,"year")) + months(7))
effort <- mutate(effort, ssn=factor(ifelse(month(Date) %in% ssn,"S","NS")),
                 occ.yy = cut(Date, breaks=occasions.yy))
sight <- mutate(sight, ssn=factor(ifelse(month(Date) %in% ssn,"S","NS")),
                occ.yy = cut(Date, breaks=occasions.yy))
indiv <- mutate(indiv, ssn=factor(ifelse(month(Date) %in% ssn,"S","NS")),
                occ.yy = cut(Date, breaks=occasions.yy))

## summarize effort types by day, by occasion
eff.yy <- effort %>% group_by(occ.yy, B, ssn) %>% summarize(Minutes=sum(Minutes,na.rm=TRUE)) %>% 
  ungroup() 
effdays <- effort %>% group_by(occ.yy) %>% summarize(ndays=n_distinct(Date))
eff.yy <- eff.yy %>% mutate(ssnB=paste(ssn, B, sep="_")) %>% 
  select(-one_of(c("ssn","B"))) %>% spread(ssnB, Minutes,fill=0) %>% 
  mutate(Min.S=rowSums(.[starts_with("S_",vars=names(.))],na.rm=T), Min.NS=rowSums(.[starts_with("NS_",vars=names(.))],na.rm=T),
         Min=Min.S+Min.NS,
         eff.02.NS=NS_0+NS_1+NS_2, eff.02.S=S_0+S_1+S_2, eff.02=eff.02.NS+eff.02.S,
         eff.03.NS=eff.02.NS+NS_3, eff.03.S=eff.02.S+S_3, eff.03=eff.03.NS+eff.03.S,
         eff.3=NS_3+S_3) %>% 
         left_join(effdays)
rm(effdays)

# NOTE: sight$lon[36] suspicious; use efflon instead if exists
#sight$lon[36] <- sight$efflon[36]
# seems unimportant even for capture distribution map


################ SELECT ZIPHIUS DATA FOR ANALYSIS #################

## index subsets of records (age class, L/R photo ID, photo quality)
indiv %<>% mutate(RQ2=qR<=2 & !is.na(qR), LQ2=qL<=2 & !is.na(qL), 
                  R2=mR>=2 & !is.na(mR), L2=mL>=2 & !is.na(mL),
                  R1=!is.na(mR), L1=!is.na(mL),
                  Q30 = (!is.na(q3R)&!q3R) | (!is.na(q3L)&!q3L), 
                  Q30L=!is.na(q3L)&!q3L, Q30R=!is.na(q3R)&!q3R,
                  nC=AgeClass!="Calf") 
## index capture occasion(s) for mR1 when RQ2&LQ2&L2 IFF at another point in CH have RQ2&R2&LQ2&L2
#indiv %<>% group_by(ID) %>% 
#  mutate(mL2R2=any(RQ2 & LQ2 & R2 & L2), mL2R1 = mL2R2 & RQ2 & LQ2 & L2 & mR==1 & !is.na(mR)) %>% 
#  ungroup() %>% select(-mL2R2)
# for each ID: any Q30 photos in ch?
indiv %<>% group_by(ID) %>% summarize(Q30.i=any(Q30)) %>% right_join(indiv)
indiv %<>% group_by(ID) %>% summarize(Q30.i.l=any(Q30L)) %>% right_join(indiv)
indiv %<>% group_by(ID) %>% summarize(Q30.i.r=any(Q30R)) %>% right_join(indiv)
# identify precaptures (captures prior to first Q30 capture)
indiv %<>% filter(Q30.i) %>% group_by(ID) %>%
  mutate(cap.Q30.occ=occ.yy[which.max(Q30)], precap=as.character(occ.yy)<as.character(cap.Q30.occ)) %>% 
  ungroup() %>% right_join(indiv) %>% mutate(precap=replace_na(precap,FALSE))
# identify individuals with linked left and right sides
# classify FinDmg (which includes ID83) as "known" histories (also those id'd from both sides)
# (ID 16 - from email chain - had distinctive fin shape and same associates 3 days in a row; only seen on one occ.yy)
indiv %<>% group_by(ID) %>% 
  mutate(L2Q2=L2&LQ2, R2Q2=R2&RQ2, linkrl=any(((L2Q2|Q30L)&(R2Q2|Q30R)) | (FinDmg&(LQ2|RQ2)) | (ID==16))) %>% 
  ungroup()

# sort indiv
indiv <- arrange(indiv, ID, Date, Sighting)

## filter data
# by photo quality only
indiv.q <- indiv %>% filter(RQ2|LQ2)
# by photo quality, distinctiveness, and age class
indiv.sel <- indiv.q %>% filter(((RQ2&R2)|(LQ2&L2)) & nC)

## data selection for analysis
# for multimark (without and with precaptures)
indiv.multi.npc <- indiv.sel %>% filter(Q30.i & !precap)
indiv.multi.pc <- indiv.sel %>% filter(Q30.i)
# Pradel left
indiv.Prad.L <- indiv.q %>% filter(nC & LQ2 & L2 & (Q30.i.l | (Q30.i.r & linkrl)))
# Pradel right
indiv.Prad.R <- indiv.q %>% filter(nC & RQ2 & R2 & (Q30.i.r | (Q30.i.l & linkrl)))

# create effort covariate that is standardized to season and Beaufort Sea State based on records with selected photo quality
indivall.q <- left_join(indiv.q, sight) %>% mutate(mm=month(Date)) %>% arrange(Date)
df.eff.cond <- effort %>% group_by(ssn,Q,B,V,S) %>% summarise(Minutes=sum(Minutes,na.rm=T)) %>% ungroup()
indivall.cond <- indivall.q %>% group_by(ssn,Q,B,V,S) %>% summarise(caps=n()) %>% ungroup()
df.cond.caps <- full_join(df.eff.cond,indivall.cond) %>% mutate(caps=replace_na(caps,0))
# beaufort and season 
df.bssn.caps  <- df.cond.caps %>% group_by(B,ssn) %>% 
  summarise(hrs=sum(Minutes,na.rm=T)/60, caps=sum(caps), cph=caps/hrs) %>% select(-hrs)
eff.bssn <- effort %>% group_by(occ.yy,B,ssn) %>% summarize(hrs=sum(Minutes,na.rm=T)/60) %>% 
  right_join(df.bssn.caps) %>% mutate(caps=hrs*cph) %>% 
  group_by(occ.yy) %>% summarize(caps=sum(caps)) %>% mutate(eff.bssn=scale(caps)) %>% 
  select(occ.yy, eff.bssn)
eff.yy <- left_join(eff.yy, eff.bssn)
rm(df.eff.cond, indivall.cond, df.cond.caps, df.bssn.caps, eff.bssn)

#table(table(indiv.Prad.L %>% group_by(ID, occ.yy) %>% summarize(cap=T) %>% select(cap, ID)))
#1  2   3  4  6 
#64 15  7  2  1 
#table(table(indiv.Prad.R %>% group_by(ID, occ.yy) %>% summarize(cap=T) %>% select(cap, ID)))
#1  2   3  4  5 
#74 22  6  2  1

## proportion population included in analysis (% sufficiently marked in photos with sufficient quality)
# can calculate at day level (use sighting level given NA IDs?) since reliably identifiable in this time frame (Wilson et al 1999)
# (no records with min Q2 photo that don't have mark grade - still holds?)
# extension of Wilson to multi-annual-level ability to re-id any animal based on admission-quality photo
# so just account for m1 animals that were left out of analysis but were admitted to catalog: 
#    once with, once without calves
# Inclusion in analysis with multimark:
# All indiv need to have same initial p regardless of m, so Q30 on one side sufficient even if m1 on that side. 
# Exclusion from analysis but same quality to be accounted for:
# Q30.i & m==1 for any sides with Q2 or better photos
# correction has to be done at annual occasion level to account for ontogenetic changes in m and AgeClass (found no stage shifts within occasion)
# one animal (149) in catalog & ! analysis, but transforms from calf to juv in later yrs
## for each ID and occasion, choose photo with highest mark score at accepted photo Q (i.e., Q30.i & (LQ2|RQ2))
# no cases where animal dips "out" of analysis population to m1 after initial inclusion
# include precaptures because make no difference to "meaningful"real" parameter outcomes for CJS, Closed
indiv.q.yy <- indiv.q %>% filter(Q30.i) %>% group_by(ID, occ.yy, AgeClass) %>% 
  summarize(n=n(), mmax = max(c(mR[!is.na(mR)&qR<=2], mL[!is.na(mL)&qL<=2])))
cf.yy <- indiv.q.yy %>% group_by(occ.yy) %>% 
  summarize(nm1woC=sum(mmax==1 & AgeClass!="Calf"), nm1plusC=sum(mmax==1 | AgeClass=="Calf"), 
            nincl=sum(mmax>1 & AgeClass!="Calf"), n=n())
cf.nc <- (1/rbeta(10000, cf.yy$nincl[1]+0.5, cf.yy$nm1woC[1]+0.5) + 1/rbeta(10000, cf.yy$nincl[2]+0.5, cf.yy$nm1woC[2]+0.5) + 
            1/rbeta(10000, cf.yy$nincl[3]+0.5, cf.yy$nm1woC[3]+0.5) + 1/rbeta(10000, cf.yy$nincl[4]+0.5, cf.yy$nm1woC[4]+0.5) + 
            1/rbeta(10000, cf.yy$nincl[5]+0.5, cf.yy$nm1woC[5]+0.5) + 1/rbeta(10000, cf.yy$nincl[6]+0.5, cf.yy$nm1woC[6]+0.5) + 
            1/rbeta(10000, cf.yy$nincl[7]+0.5, cf.yy$nm1woC[7]+0.5) + 1/rbeta(10000, cf.yy$nincl[8]+0.5, cf.yy$nm1woC[8]+0.5) + 
            1/rbeta(10000, cf.yy$nincl[9]+0.5, cf.yy$nm1woC[9]+0.5) + 1/rbeta(10000, cf.yy$nincl[10]+0.5, cf.yy$nm1woC[10]+0.5) + 
            1/rbeta(10000, cf.yy$nincl[11]+0.5, cf.yy$nm1woC[11]+0.5))/11
# sd(cf.nc)/mean(cf.nc)   # 1.9% CV, can ignore
cf.nc <- mean(cf.nc)
cf.c <- (1/rbeta(10000, cf.yy$nincl[1]+0.5, cf.yy$nm1plusC[1]+0.5) + 1/rbeta(10000, cf.yy$nincl[2]+0.5, cf.yy$nm1plusC[2]+0.5) + 
           1/rbeta(10000, cf.yy$nincl[3]+0.5, cf.yy$nm1plusC[3]+0.5) + 1/rbeta(10000, cf.yy$nincl[4]+0.5, cf.yy$nm1plusC[4]+0.5) + 
           1/rbeta(10000, cf.yy$nincl[5]+0.5, cf.yy$nm1plusC[5]+0.5) + 1/rbeta(10000, cf.yy$nincl[6]+0.5, cf.yy$nm1plusC[6]+0.5) + 
           1/rbeta(10000, cf.yy$nincl[7]+0.5, cf.yy$nm1plusC[7]+0.5) + 1/rbeta(10000, cf.yy$nincl[8]+0.5, cf.yy$nm1plusC[8]+0.5) + 
           1/rbeta(10000, cf.yy$nincl[9]+0.5, cf.yy$nm1plusC[9]+0.5) + 1/rbeta(10000, cf.yy$nincl[10]+0.5, cf.yy$nm1plusC[10]+0.5) + 
           1/rbeta(10000, cf.yy$nincl[11]+0.5, cf.yy$nm1plusC[11]+0.5))/11
# sd(cf.c)/mean(cf.c)   # 2.7% CV, can ignore
cf.c <- mean(cf.c)
#save(cf.nc, cf.c, file="cf.rdata")
rm(indiv.q.yy, cf.yy)
# alternatively, using delta method, also an inconsequential difference in uncertainty (see e.g. Sutaria and Marsh 2011),
# so just multiply distributions by multiplier


################ CREATE CAPTURE HISTORIES AND WRITE CHs AND COVARS TO FILE #################

# annual nCL
ch.yy = with(indiv.Prad.L, table(ID, occ.yy))
ch.yy[ch.yy>1] = 1
write.csv(ch.yy,file="ch.yy.Prad.L.csv",row.names=F)

# annual nCR
ch.yy = with(indiv.Prad.R, table(ID, occ.yy))
ch.yy[ch.yy>1] = 1
write.csv(ch.yy,file="ch.yy.Prad.R.csv",row.names=F)

# write annual occasion covariates to file
write.csv(eff.yy,file="occasions.cov.yy.20190315.csv",row.names=F)

rm(ch.yy)


###### PROCESSING FOR MULTIMARK #########

# create new variables for capture types
indiv.multi.npc %<>%  mutate(T1 = LQ2 & !RQ2, T2 = RQ2 & !LQ2, T4 = RQ2 & LQ2)
indiv.multi.pc %<>%  mutate(T1 = LQ2 & !RQ2, T2 = RQ2 & !LQ2, T4 = RQ2 & LQ2)

# collapse to annual occasions
indiv.multi.npc.yy <- indiv.multi.npc %>% group_by(ID, occ.yy, linkrl) %>% 
  summarize(T1=sum(T1)>0, T2=sum(T2)>0, T4=sum(T4)>0) %>% 
  mutate(captype = if_else(T4, T4*4, T1*1 + T2*2))
indiv.multi.pc.yy <- indiv.multi.pc %>% group_by(ID, occ.yy, linkrl) %>% 
  summarize(T1=sum(T1)>0, T2=sum(T2)>0, T4=sum(T4)>0) %>% 
  mutate(captype = if_else(T4, T4*4, T1*1 + T2*2))

# define vector of "known" encounter histories (same for both npc and pc)
known <- indiv.multi.npc.yy %>% group_by(ID) %>% summarize(kn=any(captype==4 | linkrl))

# create capture history
indiv.multi.npc.yy %<>%  arrange(occ.yy,ID)
indiv.multi.pc.yy %<>%  arrange(occ.yy,ID)
CH.multi.npc <- spread(indiv.multi.npc.yy[c("ID","occ.yy","captype")], key=occ.yy, value=captype, fill=0, drop=FALSE) %>% 
  ungroup() %>% select(-ID)
CH.multi.pc <- spread(indiv.multi.pc.yy[c("ID","occ.yy","captype")], key=occ.yy, value=captype, fill=0, drop=FALSE) %>% 
  ungroup() %>% select(-ID)
write.csv(CH.multi.npc, file="ch.yy.multi.npc.20190315.csv", row.names=F)
write.csv(CH.multi.pc, file="ch.yy.multi.pc.20190315.csv", row.names=F)
write.csv(known, file="multi.known.20190315.csv", row.names=F)

rm(indiv.multi.npc.yy, indiv.multi.pc.yy, CH.multi.npc, CH.multi.pc, known)

