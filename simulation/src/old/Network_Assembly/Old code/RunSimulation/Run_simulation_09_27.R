##  begin simulation

prms <- list(ntree=1,
             nrule=8,
             nmat=5,
             nloop=1,
             sp=30,
             nnull=2,
             mu=this.mu,		# from `sim_loop.R'
             lambda=this.lambda,# from `sim_loop.R'
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
 #sim.links <- mclapply(1:prms$ntree, f, mc.preschedule=FALSE)
  sim.links <- lapply(1:prms$ntree, f)
  
  res <- list(prms=prms, sim.links=sim.links)
  save(res, file=sprintf("10_13_2012/%s.RData", paste("mu",this.mu,"la",this.lambda,"rep",i,sep="_")))
}
