sim.phylo <- function(nspecies) {
  tree.1 <- rcoal(n= nspecies)
  tree.2 <- rcoal(n= nspecies)
  return(list(tree.1, tree.2))
}