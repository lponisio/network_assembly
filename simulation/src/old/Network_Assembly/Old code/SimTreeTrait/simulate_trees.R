
library(TreeSim)

sim.phylo <- function(nsim, mu, lambda, age, nspecies){
	
	tree.1 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	tree.2 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	
	return(list(tree.1, tree.2))
}

