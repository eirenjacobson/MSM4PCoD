
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(coda)


id <- "D50_LCP_Ideal1_2023-08-07"
nsim <- 10

folder <- paste0("./Figures/", id)
load(paste0("./Results/ProcResults_", id, ".RData"))

sdf <- results.out$sdf
ltdf <- results.out$ltdf
pamdf <- results.out$pamdf
rdf <- results.out$rdf
Ndf <- results.out$Ndf
bdf <- results.out$bdf
pardf <- results.out$pardf

# now plot each iteration separately
for (i in 1:nsim){
  rtest <- filter(rdf, Iter == i)
  stest <- filter(sdf, Iter == i)
  ltest <- filter(ltdf, Iter == i)
  ptest <- filter(pamdf, Iter == i)
  btest <- filter(bdf, Iter == i)

# Plot at the level of the metapopulation (A + B)
  metaplot <- ggplot()+
                geom_ribbon(data = rtest, aes(x=Year, ymin = KLCI, ymax = KUCI), fill = "lightblue")+
                geom_line(data = rtest, aes(x=Year, y = K), color = "mediumblue", lty = "dashed")+
                geom_line(data = stest, aes(x=Year, y = Ntot, group = Iter))+
                geom_point(data = ltest, aes(x=Year, y = Nhat, group = Iter))+
                geom_errorbar(data = ltest, aes(x=Year, ymin = LCI, ymax = UCI, group = Iter))+
                geom_ribbon(data = rtest, aes(x=Year, ymin = LCI, ymax = UCI), fill = "orange2", alpha = 0.5) +
                geom_line(data = rtest, aes(x=Year, y = Ntot, group = Iter), color = "orange2") +
                ggtitle("Metapopulation") +
                theme_bw()
  
  #ggsave(metaplot, filename = paste0(folder, "/metaplot_", i, ".png"), width = 5, height = 3, units = "in")

# Plot at the level of the subpopulation (B only)

  subplot <- ggplot()+
              geom_ribbon(data = rtest, aes(x=Year, ymin = KLCI/2, ymax = KUCI/2), fill = "lightblue", alpha = 0.5)+
              geom_line(data = rtest, aes(x=Year, y = K/2), color = "mediumblue", lty = "dashed")+
              geom_line(data = btest, aes(x=Year, y = N))+
              geom_point(data = ptest, aes(x=Year, y = Nhat, group = Iter), pch = 17)+
              geom_errorbar(data = ptest, aes(x=Year, ymin = LCI, ymax = UCI, group = Iter))+
    geom_ribbon(data = rtest, aes(x=Year, ymin = LCI/2, ymax = UCI/2), fill = "orange2", alpha = 0.5) +
    geom_line(data = rtest, aes(x=Year, y = Ntot/2, group = Iter), color = "orange2", alpha = 0.5) +
              ggtitle("Subpopulation") +
              theme_bw()
  
  #ggsave(subplot, filename = paste0(folder, "/subplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  # Plot of parameters to be estimated
  
  s <- pardf %>% filter(Iter == i) %>%
    pivot_longer(cols = c(K1, K2, S2, PCap), names_to = "Parameter", values_to = "Value")
  
  parplot <- ggplot(s)+
    geom_density(aes(x=Value, color = Iter), color = "mediumpurple", fill = "mediumpurple", alpha = 0.5)+
    facet_wrap(~Parameter, scales = "free") +
    theme_bw()
  
  #ggsave(parplot, filename = paste0(folder, "/parplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  chainplot <- ggplot(s) +
    geom_path(aes(x = Sample, y= Value, color = factor(Chain)), show.legend=FALSE)+
    facet_wrap(~Parameter, scales = "free") +
    theme_bw()
  
  #ggsave(chainplot, filename = paste0(folder, "/chainplot_", i, ".png"), width = 5, height = 3, units = "in")
  
  ggsave(ggarrange(metaplot, subplot, chainplot, parplot, nrow = 2, ncol = 2), 
         filename = paste0(folder, "/", id, "_", i, ".png"), 
         width = 12, height = 8, units = "in")

} # end for i
