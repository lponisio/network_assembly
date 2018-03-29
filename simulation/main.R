#rm(list=ls())
## change to your working directory, i.e., whever the github
## repository was cloned
setwd('~/pkg/network_assembly')

tree <- 'bd'

## packages to install are listed in initialize.R. Install them before
## running this line which opens all of the librarys required
source('simulation/src/all/initialize.R')

## set to the number of cores on which to run the simulation in
## parallel

options(mc.cores=1)
prms <- base.prms()

## number of replicates of each parameter combination
nreps <- 100

## brownian moition sigma options
BM <- c("BMLow", "BMHigh")

## tree depth (extinction and speciation rates) options
tree.depth <- c("Shallow", "Deep")

## loop over the parameter combinations evaluated in the study
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
        dir.create(file.path(sim.dir, type), showWarnings = FALSE)
        mclapply(1:nrow(cases), run.sim,
                 prms= prms,
                 save.path=file.path(sim.dir, type),
                 mc.preschedule = FALSE)
        dats <- load.fix(file.path(sim.dir, type), bd=tree)
        save(dats, file=file.path(sim.dir, sprintf('%s.Rdata', type)))
    }
}


## plot (figure code not on gihub at the moment)
combins <- c("BMLowDeep", "BMHighDeep", "BMLowShallow",
             "BMHighShallow")

for(i in 1:length(combins)){
    type <- combins[i]
    load(file=file.path(sim.dir, sprintf('%s.Rdata', type)))
    ## source('simulation/src/all/plotting.R')
}
