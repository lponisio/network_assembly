run.sim <- function(prms, i, save.path) {
 
  prms <- prms.f(prms)

  all.trees <- sim.phylo(nspecies=prms$sp)
  
  f <- function(x) {
    #cat(x, '\n')
   # browser()
    master.fun(prms,
               tree1 = all.trees[[1]],
               tree2 = all.trees[[2]],
               nnul = prms$nulls,
               mean.Abund = prms$mean.abund, 
               sd.Abund = prms$sd.abund,
               w.evo = prms$w.evo)
  }
  
  sim.links <- lapply(1:prms$ntree, f)

  res <- list(prms= prms, sim.links= sim.links)
  save(res, file= file.path(save.path, sprintf('%d.RData', i)))
}
