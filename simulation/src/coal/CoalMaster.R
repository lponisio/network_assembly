run.sim <- function(prms, i, save.path) {
  prms <- prms.f(prms)
  prms$range.size <- cases[i, 'range.size']
  ## prms$sp <- cases[i, 'sp']
  all.trees <- sim.phylo(nspecies=prms$sp)
  f <- function(x) {
    master.fun(prms,
               tree1 = all.trees[[1]],
               tree2 = all.trees[[2]],
               nnul = prms$nulls)
  }
  sim.links <- lapply(1:prms$ntree, f)
  res <- list(prms= prms, sim.links= sim.links)
  if( ! file.exists( save.path ) ) {
    dir.create( save.path, recursive=TRUE )
  }
  save(res, file= file.path(save.path, sprintf('%d.RData', 22036+i)))
}
