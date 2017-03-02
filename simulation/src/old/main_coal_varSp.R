rm(list=ls())
setwd("~/Dropbox/network_assembly")
source('simulation/src/all/CalcMetrics_bascompte.R')
source('simulation/src/all/load_fix.R')
source('simulation/src/all/SimTreeTrait.R')

source('simulation/src/coal/all/SimPhyloCoal.R')
source('simulation/src/coal/varNsp_master.R')
source('simulation/src/coal/all/initialize.R')

source('figures/src/mod_by_nodf/mod_by_nodf_sd.R')
source('figures/src/mod_by_nodf/mod_by_nodf_bar3.R')
source('figures/src/images/image_plot_abund_ecotype.R')
source('figures/src/corInt/corint_sd.R')
source('figures/src/singleMetrics/single_metric.R')
source('figures/src/singleMetrics/mets3.R')
source('figures/src/mods/mods.R')
source('figures/src/misc.R')

options(mc.cores= 8)
prms <- base.prms()

cases <- expand.grid(sp=seq(from=15, to=105, length=7),
                     rep=1:90)

mclapply(1:nrow(cases), run.sim,
         prms= prms,
         save.path='simulation/saved/coal/nsp',
         mc.preschedule= FALSE)

dats <- load.fix("simulation/saved/coal/nsp")

## nsp by metrics
nsp <- round(((dats[[1]][,'npol'] + dats[[1]][,'nplant'])/5))*5
dats[[1]] <- cbind(dats[[1]], nsp)

path <- "figures/coal/SingleMetrics"
mets3.plot(dats, "zNODF", "zmodD", "nsp",
           "Number of Species", path, subset=FALSE,
           column="range.size")


path <- "~/Dropbox/network_assembly/figures/coal/progression/nsp"

progress.plot(dats, "NODF", "modularityD", path)
progress.plot(dats, "zNODF", "zmodD", path)

progress.plot.cor(dats, "corPol", "corPlant", path, at=0.5)

nsp.plot(dats, "zNODF", "zmodD", "nsp", "Nestedness")

path <- "figures/coal/NODFbyMod/nsp"

mod.by.nodf.plot(simres = dats, path=path, metric1= "NODF",
                 metric2="modularityD")
mod.by.nodf.plot(simres = dats, path=path, metric1= "zNODF",
                 metric2="zmodD")

corInt.plot(simres = dats, metric1="corPol", metric2="corPlant",
            path=path)




