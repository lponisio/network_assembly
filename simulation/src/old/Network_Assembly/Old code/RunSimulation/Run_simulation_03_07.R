##  begin simulation
rm(list=ls())
setwd("~/Documents/Network_Assembly")
source("Functions/everything8.R")
source("Functions/Final_net_null.R")
library("multicore")

options(cores=24)

prms <- list(ntree=100,
             nrule=8,
             nmat=5,
             nloop=10,
             sp=30,
             mu=0.05,
             lambda=0.05,
             age=100)

for(i in 1:prms$nloop) {
  all.trees <- sim.phylo(nsim=prms$ntree,
                         mu=prms$mu,
                         lambda=prms$lambda,
                         age=prms$age,
                         nspecies=prms$sp)

  f <- function(x) {
    cat(x, "\n")
    master.fun(tree1=all.trees[[1]][[x]],
               tree2=all.trees[[2]][[x]])
  }
  sim.links <- mclapply(1:prms$ntree, f, mc.preschedule=FALSE)
  
  res <- list(prms=prms, sim.links=sim.links)
  save(res, file=sprintf("09_05_2012/%d.RData", i))	
}
