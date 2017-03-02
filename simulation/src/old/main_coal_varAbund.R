rm(list=ls())
setwd("~/Dropbox/network_assembly")
source('simulation/src/all/CalcMetrics_bascompte.R')
source('simulation/src/all/SimTreeTrait_abundSamp.R')
source('simulation/src/all/load_fix.R')

source('simulation/src/coal/all/SimPhyloCoal.R')
source('simulation/src/coal/all/initialize_abundSamp.R')
source('simulation/src/coal/varAbund_master.R')

source('figures/src/mod_by_nodf/mod_by_nodf_sd.R')
source('figures/src/images/image_plot_abund_ecotype.R')
source('figures/src/corInt/corint_sd.R')
source('figures/src/singleMetrics/single_metric.R')
source('figures/src/mods/mods.R')

source('figures/src/progression/prog_cor_multip.R')
source('figures/src/progression/prog_met_multip.R')

options(cores= 8)
prms <- base.prms()

cases <- expand.grid(mean.abund= log(seq(exp(0),exp(1), length=10)),
                     sd.abund = log(seq(exp(0),exp(1), length=10)),
                     rep=1:90)

## mclapply(1:nrow(cases), run.sim,
##          prms= prms,
##          save.path='simulation/saved/coal/bascompte',
##          mc.preschedule= FALSE)

lapply(1:nrow(cases), run.sim,
         prms= prms,
         save.path='simulation/saved/coal/misc')

dats <- load.fix("simulation/saved/coal/bascompte")

save(dats, file='simulation/saved/loaded/abundData.Rdata')

## progression plots
path <- "figures/coal/progression"

progress.plot(dats, "NODF", "modularityR", path)
progress.plot(dats, "zNODF", "zmodR", path)

progress.plot(dats, "NODF", "modularityG", path)
progress.plot(dats,  "zNODF", "zmodG", path)

progress.plot(dats, "NODF", "modularityD", path)
progress.plot(dats, "zNODF", "zmodD", path)

progress.plot.cor(dats, "corPol", "corPlant", path)

## nestedness by modularity 
## path <- "figures/coal/NODFbyMod"

path <- "figures/presentation"

mod.by.nodf.plot(simres = dats, path=path, metric1= "NODF",
                 metric2="modularityD")
mod.by.nodf.plot(simres = dats, path=path, metric1= "zNODF",
                 metric2="zmodD")

corInt.plot(simres = dats, metric1="corPol", metric2="corPlant",
            path=path)


## modularity metrics
path <- "figures/coal/mods"

mods.plot(simres = dats, metric1="modularityD", metric2="modularityR",
          metric3="modularityG", path=path)

mods.plot(simres = dats, metric1="zmodD", metric2="zmodR",
          metric3="zmodG", path=path)

## image plots
image.plot(dats, "zNODF", "eco")
image.plot(dats, "zNODF", "eco2")


