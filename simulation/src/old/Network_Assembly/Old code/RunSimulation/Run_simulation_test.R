##  begin simulation
rm(list=ls())
setwd("~/Documents/Network_Assembly")
source("Functions/everything8.R")
source("Functions/Final_net_null_09_07.R")

set.seed(1)
prms <- list(ntree=5,
             nrule=8,
             nmat=5,
             nloop=1,
             sp=30,
             mu=0.01,
             lambda=0.01,
             age=10^6)

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
  sim.links <- lapply(1:prms$ntree, f)
  
  res <- list(prms=prms, sim.links=sim.links)
  save(res, file=sprintf("09_12_2012/%d.RData", i))
}
