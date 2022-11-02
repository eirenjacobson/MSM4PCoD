############# plot capture locations and SOAR #################

# load capture data for records included in (multimark) analysis
library(dplyr)
source("1_ziphius.process.pub.r")
indivall.multi.pc <- left_join(indiv.multi.pc, sight) %>% select(lat, lon)

#install.packages("digest")
#library(devtools)
#install_github("tidyverse/ggplot2")
#devtools::install_github('jmlondon/ptolemy')
#install_gshhg()
library(ptolemy)   # calcur
library(ggplot2)
library(dplyr)
library(sf)   #st_transform, as_Spatial
library(readr)   # read_csv
library(marmap)   # autoplot.bathy
library(maps)
library(mapdata)

# set coordinates for immediate and larger-scale bounding boxes
map.sm <- data.frame(lat=c(32.5, 33.3), lon=c(-119.35, -118.1))
map.lg <- data.frame(lat=c(32.25, 35.4), lon=c(-121, -116.5))
aspect <- cos((mean(map.sm$lat)/90)*(pi/2))
aspectlg <- cos((mean(map.lg$lat)/90)*(pi/2))

# bathymetry 
snbathy <- read_csv("../../ERDDAP_bathy_SNB.csv", col_types = "nnn", skip=2, col_names=c("lat","lon","z")) %>% 
  filter(lat>map.sm$lat[1], lat<map.sm$lat[2], lon<map.sm$lon[2], lon>map.sm$lon[1]) %>% select(lon,lat,z) %>% 
  mutate(z = z*as.numeric(z<0))
# turn into rectangular format for compatibility with marmap
snbathy <- snbathy %>% arrange(lat, lon)
z <- matrix(snbathy$z, ncol=length(unique(snbathy$lat)), 
            dimnames=list(unique(snbathy$lon), unique(snbathy$lat)))
class(z) <- "bathy"

# SOAR polygon
soarsf <- st_read("SOAR.KML")
soar_ll <- data.frame(matrix(unlist(soarsf[1,]$geometry),ncol=3)) %>% 
  rename(lon=X1, lat=X2) %>% select(lon,lat)

# get coastline data, convert to latlon
calcur_base_h <- ptolemy::calcur(resolution="h", simplify=F) %>% 
  sf::st_transform(4326) %>% filter(!is.na(st_dimension(geometry)))
coast <- fortify(as_Spatial(calcur_base_h)) %>% 
  filter(lat>map.sm$lat[1], lat<map.sm$lat[2], long<map.sm$lon[2], long>map.sm$lon[1])

mainmap <- autoplot.bathy(z, geom = "tile") +
  scale_fill_gradient(low=gray(2/7), high=gray(6.5/7)) + 
  theme(legend.position="none", plot.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text=element_text(size=8)) + 
  labs(x=NULL, y=NULL) + 
  scale_x_continuous(labels = c("","119.1°W","118.8°W","118.5°W","118.2°W"), expand = c(0,0)) + 
  scale_y_continuous(labels = label_number(accuracy = 0.1, suffix = "°N"), expand = c(0,0)) +
  theme(panel.border = element_rect(linetype = 1, fill = NA)) +
  geom_contour(data=snbathy, aes(x=lon,y=lat,z=z), breaks = c(-800, -1200, -1600), col="grey24", size = 0.15) +
  geom_label(data = data.frame(x=c(-118.32, -118.23), y=c(32.75, 32.65), 
                               label = c("800 m", "1600 m")), 
             aes(x = x, y = y, label = label), size = 2.25, label.padding = unit(0.05, "lines"), 
             label.size = 0, fill = "white") +
  geom_polygon(data = coast, aes(x=long, y=lat, group=group), fill = "white", col="black", size=0.2) + 
  annotate("text", label="San Clemente Isl", x=-118.489, y=32.908, angle=309.5, size = 2.25) +
  geom_polygon(data=soar_ll, aes(x=lon, y=lat), color="white", size=1, fill = NA) +
  geom_point(data=indivall.multi.pc, aes(x=lon, y=lat), shape=24, size=1.5, fill="white")

# inset map
usa <- map_data("worldHires", region=c("USA","Mexico","California"))
insetmap <- ggplotGrob( ggplot() + 
                          geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
                                       fill = "grey80", colour="black", size=0.2) + 
                          theme(panel.background = element_rect(fill = "white"), 
                                axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),
                                axis.title.x=element_blank(), axis.title.y=element_blank(),
                                axis.text = element_blank(), axis.ticks.length = unit(0, "mm"),
                                panel.border = element_rect(linetype = 1, fill = NA), 
                                plot.margin = margin(0,0,0,0)) +
                          annotate("rect", xmin=map.sm$lon[1], xmax=map.sm$lon[2], 
                                   ymin=map.sm$lat[1], ymax=map.sm$lat[2], alpha=0, color="black", size=0.5) +
                          annotate("text", x=-118.8, y=34.9, label="CALIFORNIA", hjust=0, angle=322, size=2) +
                          coord_fixed(1/aspectlg, xlim = map.lg$lon, ylim=map.lg$lat)
)

# create figure with inset map
g3 <- mainmap +
  annotation_custom(grob = insetmap, xmin = -118.47, xmax = map.sm$lon[2]-0.002,
                    ymin=33.04, ymax = map.sm$lat[2])

ggsave(plot=g3, filename="Fig1_Map_captures.tiff", 
       width=110, height = 80, units="mm", dpi = 1200, compression = "lzw")


############# plot total effort, BSS 0-2 effort by season and year #################

library(lubridate)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

eff.yy <- read_csv("occasions.cov.yy.20190315.csv")

eff.yy %<>% select(occ.yy, Min.S, Min.NS) %>% 
  pivot_longer(cols = starts_with("Min"), names_to = "ssn", names_prefix = "Min.", 
               values_to = "Min", values_drop_na = FALSE) %>% 
  rowwise() %>% 
  mutate(occ.yy = year(occ.yy), hr = Min/60, 
         ssn = switch(ssn, S = "Summer", NS = "Nonsummer")) %>% 
  ungroup() %>% 
  filter(Min>0)

p <- ggbarplot(eff.yy, x = "occ.yy", y = "hr", 
          fill = "ssn", color = "ssn",
          palette = c("black", "grey60"),
          xlab="Occasion", ylab="Search effort (hours)",
          legend.title = "",
          x.text.angle = 300)
pcust <- ggpar(p, font.tickslab = 8, font.legend = 8, font.x = 10, font.y = 10) +
  theme(legend.key.size = unit(0.7, "lines"), 
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "line"), hjust = -0.1),
        axis.title.x = element_text(margin = margin(t = 0.5, unit = "line")),
        legend.position = c(0.5, 0.975), legend.direction = "horizontal",
        plot.margin = margin(1, 0.2, 0, 0.2))

ggsave(plot=pcust, filename="Fig3_effort.tiff", 
       width=80, height = 55, units="mm", dpi = 1200, compression = "lzw")


################# plot caps/hr vs beaufort, season ####################

library(ggplot2)
library(magrittr)

source("1_ziphius.process.pub.r")

# create effort covariate that is standardized to season and Beaufort Sea State based on records with selected photo quality
indivall.q <- left_join(indiv.q, sight) %>% mutate(mm=month(Date)) %>% arrange(Date)
df.eff.cond <- effort %>% group_by(ssn,Q,B,V,S) %>% summarise(Minutes=sum(Minutes,na.rm=T)) %>% ungroup()
indivall.cond <- indivall.q %>% group_by(ssn,Q,B,V,S) %>% summarise(caps=n()) %>% ungroup()
df.cond.caps <- full_join(df.eff.cond,indivall.cond) %>% mutate(caps=replace_na(caps,0))
df.bssn.caps  <- df.cond.caps %>% group_by(B,ssn) %>% 
  summarise(hrs=sum(Minutes,na.rm=T)/60, caps=sum(caps), cph=caps/hrs) %>% rename(Season=ssn) %>% 
  mutate(Season = if_else(Season == "S", "Summer", "Nonsummer"))

g1 <- ggplot(df.bssn.caps, aes(B, cph)) + 
  geom_point(aes(colour=Season, shape=Season, size=hrs), size=2) + 
  scale_colour_manual(values=c("black","gray60"), name = "") + 
  scale_shape_discrete(name = "") + 
  xlab("Beaufort Sea State") + ylab("Captures per Hour") +
  theme_classic() + theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 8), 
                          legend.position = c(0.8, 0.95), plot.margin = margin(1, 0.2, 0, 0.2))

ggsave(plot=g1, filename="Fig4_capsperhr_bssn.tiff", 
       width=80, height=60, units="mm", dpi = 1200, compression = "lzw")


###################### plot discovery curves ########################

### set up on annual basis because presence in SOAR autocorrelated within sampling periods based on tag data
# set up annual ch, then growing ch for each year, and get caps vs indiv
library(readr)
ch.yy.l  <- as.matrix(read_csv("ch.yy.Prad.L.csv"))
ch.yy.r  <- as.matrix(read_csv("ch.yy.Prad.R.csv"))

cumu.r = data.frame(ids=integer(),ind=integer())
cumu.l = data.frame(ids=integer(),ind=integer())
for (y in 1:dim(ch.yy.r)[2]) {
  temp.r = matrix(ch.yy.r[,1:y],ncol=y)
  temp.r = matrix(temp.r[rowSums(temp.r)>0,],ncol=y)
  cumu.r = rbind(cumu.r,data.frame(ids=sum(temp.r),ind=dim(temp.r)[1]))
  temp.l = matrix(ch.yy.l[,1:y],ncol=y)
  temp.l = matrix(temp.l[rowSums(temp.l)>0,],ncol=y)
  cumu.l = rbind(cumu.l,data.frame(ids=sum(temp.l),ind=dim(temp.l)[1]))
}

opar = par()

tiff(filename="Fig5_discovery.tiff", width=100, height=75, units="mm", pointsize=10, res=1200, compression="lzw")

par(mgp=c(2,0.6,0), mar = c(3,3,0.1,0.1)) 
plot(c(0,cumu.r$ids),c(0,cumu.r$ind),type="l",col=1,xlab="Cumulative identifications",ylab="Cumulative individuals",
     xlim=c(0,155),ylim=c(0,155),lwd=2,xaxs="i",yaxs="i",cex.lab=1, cex.axis=0.8)
lines(c(0,cumu.l$ids),c(0,cumu.l$ind),col="gray60",lwd=2)
points(cumu.r$ids,cumu.r$ind,col=1,cex=1.25)
points(cumu.l$ids,cumu.l$ind,col="gray60",pch=2,cex=1.25)
abline(0,1,lty=2)

dev.off()

par = opar

#################### plot posteriors ##################

# plot posteriors for phi, N, lambda, and b

library(ggplot2)
library(coda)
library(gridExtra)

attach("multimark.cjs.final.rdata")
phi = data.frame(x=as.matrix(multi.effbssniRE.pc.long.p[,11])[,1],
                 b=as.matrix(multi.effbssniRE.pc.long$mcmc[,2])[,1])
detach(2)

attach("multimark.closed.20190315.rdata")
attach("cf.rdata")   # total-population correction factor from ziphius.selectdata.r
N = data.frame(x=cf.c*as.matrix(multi.effiRE.sub$mcmc[,"N"])[,1])
detach(2)
detach(2)

attach("PradelLambda.20190315.rdata")
lam = data.frame(x=as.matrix(draws.PL.pefftRE.phifix.R[,"mean.rho"])[,1])
detach(2)

g1 <- ggplot(data=phi, aes(x)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab("Apparent annual survival") + ylab("Density") +
  scale_x_continuous(limits = c(0.8, 1)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_text(aes(0.77, -1, label = "Probability density"))

g2 <- ggplot(data=N, aes(x)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab("Population size") + ylab("Density") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 450)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

g3 <- ggplot(data=phi, aes(b)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab("Effort coefficient") + geom_vline(xintercept=0, linetype=2) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.2,0.6)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

g4 <- ggplot(data=lam, aes(x)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab("Population growth rate") + 
  scale_x_continuous(limits = c(0.85,1.15)) + 
  geom_vline(xintercept=1, linetype=2) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave(plot=marrangeGrob(grobs = list(g1, g2, g3, g4), nrow = 2, ncol=2, 
                         left = textGrob("Probability density", gp=gpar(fontsize=10), rot = 90), top = NULL), 
       filename="Fig6_posteriors.tiff", 
       width=180, height=100, units="mm", compression = "lzw")


################# plot inference error simulation results #########################

source("9_sims.Bayes.Pradel.postprocess.pub.r")

# plot inference error curves for current effort
opar = par()

tiff(filename="Fig7_inferror_current.tiff", width=100, height=67, units="mm", pointsize=10, res=1200, compression="lzw")
par(mai=c(0.5,0.5,0.02,0), mgp=c(2,0.6,0))

plot(lam1.error1$a, lam1.error1$`error I`, type="l",lty=1, lwd=2,
     xlim=c(0.475,0.975), ylim=c(0,0.9), xaxs="i", yaxs="i", cex.axis = 0.8,
     xlab=expression(paste(plain('Evidence threshold, 1 - Quasi-'),symbol(alpha))), ylab="Error rate")
lines(lam.dec.error2$a, lam.dec.error2$`1`, lty=2, lwd=2, col="gray70")
lines(lam.dec.error2$a, lam.dec.error2$`2`, lty=4, lwd=2, col="gray70")
lines(lam.dec.error2$a, lam.dec.error2$`15`, lty=3, lwd=2, col="gray70")
points(lam1.error1$a, lam1.error1$`error I`, cex=1.5)
points(lam.dec.error2$a, lam.dec.error2$`1`, col="gray70", pch=2, cex=1.5)
points(lam.dec.error2$a, lam.dec.error2$`2`, col="gray70", pch=5, cex=1.5)
points(lam.dec.error2$a, lam.dec.error2$`15`, col="gray70", pch=22, cex=1.5)
legend("topleft", bty = "n", y.intersp = 0.9, 
       c(expression(paste(symbol(lambda),plain(' = 1'))), expression(paste(symbol(lambda),plain(' = 0.966'))),
         expression(paste(symbol(lambda),plain(' = 0.933'))),expression(paste(symbol(lambda),plain(' = 0.87')))), 
       title = "Population growth rate", lty=c(1,2,4,3), pch=c(1,2,5,22), col=c(1,rep("gray70",3)), cex=0.9)

dev.off()

par(opar)


# plot inference error for additional effort scenarios
opar = par()

tiff(filename="Fig8_inferror_future.tiff", width=100, height=67, units="mm", pointsize=10, res=1200, compression="lzw")
par(mai=c(0.5,0.5,0.02,0), mgp=c(2,0.6,0))

ny = c(11,16,21)
# more years
plot(ny, sims.stats[c(2,9,5),"dec.80"], ylim=c(0.4,1), xaxt="n", xlab="Years of effort",
     ylab="Probability of detecting decline", cex=1.5, pch=16, cex.axis = 0.8)
axis(1, ny, labels=ny, cex.axis = 0.8)
lines(ny, sims.stats[c(2,9,5),"dec.80"], lwd=2)
points(ny, sims.stats[c(4,10,6),"dec.80"], cex=1.5,col="gray70",pch=16)
lines(ny, sims.stats[c(4,10,6),"dec.80"],col="gray70", lwd=2)

# more years and more effort
points(ny, sims.stats[c(7,13,11),"dec.80"],cex=1.5,pch=2)
lines(ny, sims.stats[c(7,13,11),"dec.80"], lty=2, lwd=2)
points(ny, sims.stats[c(8,14,12),"dec.80"], cex=1.5,col="gray70",pch=2)
lines(ny, sims.stats[c(8,14,12),"dec.80"],col="gray70",lty=2, lwd=2)

legend("bottomright",legend=c(expression(paste(lambda," = 0.966, current effort")),
                              expression(paste(lambda," = 0.933, current effort")),
                              expression(paste(lambda," = 0.966, 2x effort")),
                              expression(paste(lambda," = 0.933, 2x effort"))),
       pch=c(16,16,2,2), lty=c(1,1,2,2), col=c("black","gray70","black","gray70"), cex=0.9, pt.cex=1,
       title = "Growth rate and effort scenario", bty = "n", y.intersp = 0.9, inset = c(0.01, -0.02))
dev.off()

par(opar)
