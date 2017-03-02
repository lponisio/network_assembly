rm(list=ls())
setwd("~/Dropbox/network_assembly")
source('simulation/src/coal/initialize_coal.R')
source('simulation/src/all/CalcMetrics.R')
source('simulation/src/coal/SimTreeTrait_coal.R')
source('simulation/src/coal/load_fix_coal.R')
source('simulation/src/coal/master_coal.R') 
source('figures/code/srcparam2_plots.R')
source('figures/code/mod_by_nodf_byRule.R')

options(cores= 10)
prms <- base.prms()
#lapply(1:100, run.sim,
#      prms=prms,
#      save.path='simulation/saved/coal')
 
mclapply(1:1000, run.sim,
         prms= prms,
         save.path='simulation/saved/coal',
         mc.preschedule= FALSE)

dats <- load.fix("simulation/saved/coal",
                 "simulation/saved/coal",
                 prms$nrule, prms$nmat, 4)

path <- "~/Dropbox/network_assembly/figures/coal"

Rparam.plot(simres = dats, metric = 'NODF', ymax=100, path=path)

Rparam.plot(simres = dats, metric = 'modularity', ymax=1, path=path)

Rparam.plot(simres = dats, metric = 'Shannon diversity', ymax=10, path=path)

Rparam.plot(simres = dats, metric = 'interaction evenness', ymax=2, path=path)

Rparam.plot(simres = dats, metric = 'H2', ymax=1, path=path)

mod.by.nodf.plot(simres = dats, path=path)



