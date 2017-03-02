
library(ape)
library(mvtnorm)

character.evo <- function(tree1, tree2){
	
	tree1.trait1 <- exp(rTraitCont(tree1, model="BM"))
	tree1.trait2 <- exp(rTraitCont(tree1, model="BM"))
	
	tree2.trait1 <- exp(rTraitCont(tree2, model="BM"))
	
	sample.trait <- sample(log(tree1.trait1), size = Ntip(tree2), replace=TRUE)
	
	tree2.trait2 <- as.vector(exp(rmvnorm(1,sample.trait,sigma=0.05*vcv(tree2),method="chol")))
	
	same.same <- list(tree1.trait1, tree1.trait1)
	same.diff <- list(tree1.trait1, tree1.trait2)		
	diff.diff <- list(tree1.trait1, tree2.trait1)
	diff.same <- list(tree1.trait1, tree2.trait2)	
	
	
	return(list(same.same = same.same, same.diff = same.diff, diff.diff = diff.diff, diff.same = diff.same))
}





