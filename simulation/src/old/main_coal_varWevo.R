rm(list=ls())
setwd("~/Dropbox/network_assembly")
source('simulation/src/coal/initialize_coal.R')
source('simulation/src/all/CalcMetrics.R')
source('simulation/src/coal/SimTreeTrait_poilog.R')
source('simulation/src/coal/load_fix_coal.R')
source('simulation/src/coal/master_coal_wevo.R') 
source('figures/src/Rparam2_plots.R')
source('figures/src/image_plot_abund.R')
source('figures/src/mod_by_nodf_quantile.R')
source('figures/src/corint.R')
source('figures/src/wevo_by_metric_eco.R')

options(cores= 15)
prms <- base.prms()

cases <- expand.grid(w.evo=exp(seq(log(0.01),log(1000),length=10)),
                    rep=1:10)

#save(cases, file='simulation/saved/coal/cases.RData')

## lapply(1:nrow(cases), run.sim,
##      prms=prms,
##      save.path='simulation/saved/coal/test')
 
mclapply(1:nrow(cases), run.sim,
         prms= prms,
         save.path='simulation/saved/coal/110513_evo_null',
         mc.preschedule= FALSE)

dats <- load.fix("simulation/saved/coal/110513_evo_null",
                 "simulation/saved/coal/110513_evo_null",
                 prms$nrule, prms$nmat, 4)


mod.by.nodf.plot(simres = dats, path=path, metric1= "NODF",
                 metric2="modularity")
mod.by.nodf.plot(simres = dats, path=path, metric1= "z.NODF", metric2="z.mod")
mod.by.nodf.plot(simres = dats, path=path, metric1= "z2.NODF", metric2="z2.mod")
corInt.plot(simres = dats, path=path)
w.evo.plot(simres = dats, metric= "NODF")
w.evo.plot(simres = dats, metric= "modularity")
w.evo.plot(simres = dats, metric= "H2")
w.evo.plot(simres = dats, metric= "z.NODF")
w.evo.plot(simres = dats, metric= "z.mod")
w.evo.plot(simres = dats, metric= "z.H2")

mod.by.nodf.plot(simres = dats, path=path, metric1= "NODF",
                 metric2="modularity")
mod.by.nodf.plot(simres = dats, path=path, metric1= "z.NODF", metric2="z.mod")
mod.by.nodf.plot(simres = dats, path=path, metric1= "z2.NODF", metric2="z2.mod")
corInt.plot(simres = dats, path=path)
w.evo.plot(simres = dats, metric= "NODF")
w.evo.plot(simres = dats, metric= "modularity")
w.evo.plot(simres = dats, metric= "H2")
w.evo.plot(simres = dats, metric= "z.NODF")
w.evo.plot(simres = dats, metric= "z.mod")
w.evo.plot(simres = dats, metric= "z.H2")
