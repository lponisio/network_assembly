rm(list=ls())
## change to your working directory
setwd('~/Dropbox/network_assembly')
tree <- 'bd'
source('simulation/src/all/initialize.R')

options(mc.cores=1)
prms <- base.prms()
prms$sp <- 15
nreps <- 100

BM <- c("BMLow", "BMHigher")
tree.depth <- c("Shallow", "Deep")

for(i in BM){
    if(i == "BMLow"){
        prms$sigma <- 0.01
    } else{
        prms$sigma <- 4
    }
    for(j in tree.depth){
        if(j == "Shallow"){
            mu <- 0.01
            lambda <- 0.05
        } else{
            mu <- 0.5
            lambda <- 10

        }
        type <- paste(i, j, sep="")
        print(type)
        cases <- expand.grid(mu= mu,
                             lambda= lambda,
                             range.size= seq(0.1, 1, length=5),
                             rep=1:nreps)

        mclapply(1:nrow(cases), run.sim,
                 prms= prms,
                 save.path=file.path(sim.dir, type),
                 mc.preschedule = FALSE)
        dats <- load.fix(file.path(sim.dir, type), bd=tree)
        save(dats, file=file.path(sim.dir, sprintf('%s.Rdata', type)))
    }
}


## plot all the different combinations
combins <- c("BMLowDeep", "BMHighDeep", "BMLowShallow",
             "BMHighShallow")

for(i in 1:length(combins)){
    type <- combins[i]
    dats <- load.fix(file.path(sim.dir, type), bd=tree)
    load(file=file.path(sim.dir, sprintf('%s.Rdata', type)))
    source('simulation/src/all/plotting.R')
}
