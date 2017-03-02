rm(list=ls())
setwd("~/Dropbox/network_assembly")
source('simulation/src/all/CalcMetrics_bascompte.R')
source('simulation/src/all/load_fix.R')
source('simulation/src/all/SimTreeTrait.R')

source('simulation/src/bd/SimPhyloBD.R')
source('simulation/src/bd/bd_master.R')
source('simulation/src/bd/initialize.R')

source('figures/src/mod_by_nodf/mod_by_nodf_sd.R')
source('figures/src/mod_by_nodf/mod_by_nodf_bar3.R')
source('figures/src/images/image_plot_abund_ecotype.R')
source('figures/src/corInt/corint_sd.R')
source('figures/src/singleMetrics/single_metric.R')
source('figures/src/singleMetrics/mets3.R')
source('figures/src/intimacy/intimacyNestModCon.R')
source('figures/src/mods/mods.R')
source('figures/src/misc.R')

source('figures/src/progression/prog_cor_multip.R')
source('figures/src/progression/prog_met_multip.R')


options(mc.cores= 3)

prms <- base.prms()

cases <- expand.grid(mu= seq(from= 0.01, to= 0.5, length.out= 10),
                     lambda= seq(from= 0.05, to= 10, length.out= 10),
                     age= seq(from= 1, to = 100, length.out= 1),
                     rep=1:8)
##save(cases, file='simulation/saved/cases.RData')

mclapply(1:nrow(cases), run.sim,
         prms= prms,
         cases=cases,
         save.path='simulation/saved/bd/varTree',
         mc.preschedule= FALSE)

dats <- load.fix("simulation/saved/bd/bascompte")


## progression
path <- "figures/bd/progression"


progress.plot(dats, "NODF", "modularityR", path)
progress.plot(dats, "zNODF", "zmodR", path)

progress.plot(dats, "NODF", "modularityG", path)
progress.plot(dats,  "zNODF", "zmodG", path)

progress.plot(dats, "NODF", "modularityD", path)
progress.plot(dats, "zNODF", "zmodD", path)

progress.plot.cor(dats, "corPol", "corPlant", path)

## mod by nodf
path <- "figures/bd/NODFbyMod"

mod.by.nodf.plot(simres = dats, path=path, metric1= "NODF",
                 metric2="modularityD")
mod.by.nodf.plot(simres = dats, path=path, metric1= "zNODF",
                 metric2="zmodD")

corInt.plot(simres = dats, metric1="corPol", metric2="corPlant",
            path=path)


## image plots
image.plot(simres = dats, metric = 'NODF', 'evo')
image.plot(simres = dats, metric = 'NODF', 'eco')
image.plot(simres = dats, metric = 'NODF', 'eco2')

image.plot(simres = dats, metric = 'modularity', 'evo')
image.plot(simres = dats, metric = 'modularity', 'eco')
image.plot(simres = dats, metric = 'modularity', 'eco2')

image.plot(simres = dats, metric = 'H2', 'evo')
image.plot(simres = dats, metric = 'H2', 'ecos')

image.plot(simres = dats, metric = 'z.NODF', 'evo')
image.plot(simres = dats, metric = 'z.NODF', 'ecos')

image.plot(simres = dats, metric = 'z.mod.d', 'evo')
image.plot(simres = dats, metric = 'z.mod.d', 'ecos')

image.plot(simres = dats, metric = 'z.H2', 'evo')
image.plot(simres = dats, metric = 'z.H2', 'ecos')



