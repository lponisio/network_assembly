##  needed libraries
library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)
library(distr)



###simulate trees

sim.phylo <- function(nsim, mu, lambda, age, nspecies){
	
	tree.1 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	tree.2 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	
	return(list(tree.1, tree.2))
}


##simulate character evolution
character.evo <- function(tree1, tree2){
	
	tree1.trait1 <- exp(rTraitCont(tree1, model="BM"))
	tree1.trait2 <- exp(rTraitCont(tree1, model="BM"))
	
	tree2.trait1 <- exp(rTraitCont(tree2, model="BM"))
	
	sample.trait <- sample(log(tree1.trait1), size = Ntip(tree2), replace=TRUE)
	
	tree2.trait2 <- as.vector(exp(rmvnorm(1,sample.trait,sigma=0.05*vcv(tree2),method="chol")))
	
	## phylogenetic distances between tips in plant and pol phylos
	tree1.dist <- cophenetic(tree1)
	tree2.dist <- cophenetic(tree2)
	
	same.same <- list(tree1.trait1, tree1.trait1, tree1.dist, tree1.dist)
	same.diff <- list(tree1.trait1, tree1.trait2, tree1.dist, tree1.dist)		
	diff.diff <- list(tree1.trait1, tree2.trait1, tree1.dist, tree2.dist)
	diff.same <- list(tree1.trait1, tree2.trait2, tree1.dist, tree2.dist)
	
	return(list(same.same = same.same, same.diff = same.diff, 
				diff.diff = diff.diff, diff.same = diff.same))
}


##takes overlap matrices and abundance matrix and calculates the overlap*abundance and sampled matrices


interact.sample <- function(abundance, trait.overlap, npol, nplant){
	
	trait <- matrix(trait.overlap, nrow=npol)
	
	x <- matrix(0, nrow=npol, ncol=nplant)
	
	if(sum(trait) != 0 & length(trait[trait != 0]) >= 4){  
	
		trait.abund <- trait*abundance
	
		fill <- round(sum(trait.abund))
	
		p.trait.abund <- trait.abund/sum(apply(trait.abund,1,sum))
	
		trait.samp <- try(matrix(rmultinom(1, fill, prob=as.vector(p.trait.abund)), nrow=nrow(p.trait.abund)), silent=TRUE)
		trait.samp100 <- try(matrix(rmultinom(1, 100, prob=as.vector(p.trait.abund)), nrow=nrow(p.trait.abund)), silent=TRUE)
		trait.samp1000 <- try(matrix(rmultinom(1, 1000, prob=as.vector(p.trait.abund)), nrow=nrow(p.trait.abund)), silent=TRUE)
		
		if(class(trait.samp) == "try-error") {
	  		return(list(x,x,x,x,x))
	  	
	  	} else{
	  	
	  		return(list(trait, trait.abund, trait.samp, trait.samp100, trait.samp1000))
		}
	
	} else{
		
		return(list(x,x,x,x,x))
	}	
}

##  takes two trait vectors (simulated from phylogenies) and computes
##  different plant-animal matrices based on different linkage rules
make.pa.mat <- function(traits) {
	
	npol <- nplant <- length(traits[[1]])
	pa.comb <- expand.grid(1:npol , 1:nplant)
	
	pol.trait <- traits[[1]] 
	plant.trait <- traits[[2]]
	
	##  continuous with trait overlap

	# wide 
	range.pol.wide <- cbind((pol.trait - exp(pol.trait)/2), (pol.trait + exp(pol.trait)/2))
	range.pol.wide[range.pol.wide < 0] <- 0.01
	
	range.plant.wide <-  cbind((plant.trait - exp(plant.trait)/2), (plant.trait + exp(plant.trait)/2))
	range.plant.wide[range.plant.wide < 0] <- 0.01

	# narrow
	range.pol.narrow <- cbind((pol.trait - exp(-pol.trait)/2), (pol.trait + exp(-pol.trait)/2))
	range.pol.narrow[range.pol.narrow < 0] <- 0.01
	
	range.plant.narrow <-  cbind((plant.trait - exp(-plant.trait)/2), (plant.trait + exp(-plant.trait)/2))
	range.plant.narrow[range.plant.narrow < 0] <- 0.01
	
	## narrow traits 
	pol.narrow.min <- range.pol.narrow[pa.comb[,1],1]
	pol.narrow.max <- range.pol.narrow[pa.comb[,1],2]
	
	plant.narrow.min <- range.plant.narrow[pa.comb[,2],1]
	plant.narrow.max <- range.plant.narrow[pa.comb[,2],2]

	## wide traits
	pol.wide.min <- range.pol.wide[pa.comb[,1],1]
	pol.wide.max <- range.pol.wide[pa.comb[,1],2]
	
	plant.wide.min <- range.plant.wide[pa.comb[,2],1]
	plant.wide.max <- range.plant.wide[pa.comb[,2],2]

	##simulate abundances and define interaction probability 

	abund.pol <- rlnorm(npol, meanlog=1, sdlog=1)
	abund.plant <- rlnorm(nplant, meanlog=1, sdlog=1)
	
	abund <- outer(abund.plant, abund.pol)
	
	##  narrow overlap
	
	lower.lim.narrow <- ifelse(pol.narrow.min > plant.narrow.min, pol.narrow.min, plant.narrow.min)
	upper.lim.narrow <- ifelse(pol.narrow.max < plant.narrow.max, pol.narrow.max, plant.narrow.max)
	
	trait.overlap.narrow <- upper.lim.narrow - lower.lim.narrow 
	trait.overlap.narrow[trait.overlap.narrow < 0] <- 0
	
	trait.narrow <- interact.sample(abund, trait.overlap.narrow, npol, nplant)
	
	##  wide overlap
	
	lower.lim.wide <- ifelse(pol.wide.min > plant.wide.min, pol.wide.min, plant.wide.min)
	upper.lim.wide <- ifelse(pol.wide.max < plant.wide.max, pol.wide.max, plant.wide.max)
	
	trait.overlap.wide <- upper.lim.wide - lower.lim.wide 
	trait.overlap.wide[trait.overlap.wide < 0] <- 0
	
	trait.wide <- interact.sample(abund, trait.overlap.wide, npol, nplant)
	##  one wide, one narrow
	
	lower.lim.comp <- ifelse(pol.wide.min > plant.narrow.min, pol.wide.min, plant.narrow.min)
	upper.lim.comp <- ifelse(pol.wide.max < plant.narrow.max, pol.wide.max, plant.narrow.max)
	
	trait.overlap.comp <- upper.lim.comp - lower.lim.comp 
	trait.overlap.comp[trait.overlap.comp < 0] <- 0
	
	trait.comp <- interact.sample(abund, trait.overlap.comp, npol, nplant)	
	
	## continuous with trait barrior
	
	plant.trait.combin <- plant.trait[pa.comb[,2]]
	pol.trait.combin <- plant.trait[pa.comb[,1]]
	
	barrior <- ifelse(pol.trait.combin > plant.trait.combin, 1, 0)
	
	trait.bar <- interact.sample(abund, barrior, npol, nplant)
	## equal probability of interating, matrix of 1s
	
	##combination of barrior and comp
	trait.bar.comp <- interact.sample(abund, barrior*trait.overlap.comp, npol, nplant)
	trait.bar.narrow <- interact.sample(abund, barrior*trait.overlap.narrow, npol, nplant)
	trait.bar.wide <- interact.sample(abund, barrior*trait.overlap.wide, npol, nplant)
	
	neutral <- interact.sample(abund, matrix(1, ncol=nplant, nrow=npol), npol, nplant)

	 interact <- c(
	 	trait.comp = trait.comp, 
	 	trait.bar = trait.bar,
	 	neutral = neutral, 
	 	trait.narrow = trait.narrow,
	 	trait.wide = trait.wide,
	 	trait.bar.comp = trait.bar.comp,
	 	trait.bar.narrow = trait.bar.narrow,
	 	trait.bar.wide = trait.bar.wide )
	 	
	
	return(interact)
}

##  uses `make.pa.mat' to produce a list of `fundamental' matrices
##  network statistics are then calculated on those matrices and returned

link.rule <- function(traits) {
	
	tree.dist <- traits[3:4]
	
	interact <- make.pa.mat(traits[1:2])
	
	##  calculate metrics for each matrix in fundamental
	out.metrics <- sapply(interact,network.metrics, N=2)
	
	##  calculate phylo signal of interactions
	out.cors <- sapply(interact, function(x) interact.phylo.signal(tree.dist,x))
	
	
	return(rbind(out.metrics,out.cors))
	
}

##calculate phylogenetic structure of interactions
##tree 1 = polinators, tree 2 = plants, int.mat is a list of interaction matrices

interact.phylo.signal <- function(dist.tree, int.mat) {
	dist.tree1 <- dist.tree[[1]]
	dist.tree2 <- dist.tree[[2]]
	
	print(int.mat)
	print(sum(int.mat) != 0 & all(is.na(int.mat) == FALSE))
	
	if(sum(int.mat) != 0 & all(is.na(int.mat) == FALSE)){
		
		cNonEmpty <- colSums(int.mat) > 0
		rNonEmpty <- rowSums(int.mat) > 0
		
		int.mat <- int.mat[rNonEmpty,cNonEmpty]
		dist.tree1 <- dist.tree1[cNonEmpty,cNonEmpty]
    	dist.tree2 <- dist.tree2[rNonEmpty, rNonEmpty]
    	
    	print(int.mat)
    	print(is.matrix(int.mat))
    	print(all(dim(int.mat) >= 10))
    	
    	if(is.matrix(int.mat)){
    		
    		if(all(dim(int.mat) >= 10)){
    		
				plant.dist <- vegdist(int.mat, method="bray", diag=TRUE, upper=TRUE) 
				pol.dist <- vegdist(t(int.mat), method="bray", diag=TRUE, upper=TRUE)
		
				cor.plant <- mantel(plant.dist, dist.tree2, permutations=FALSE, method="spearman")$statistic
				cor.pol <- mantel(pol.dist, dist.tree1, permutations=FALSE, method="spearman")$statistic
	
				return(c(cor.pol = cor.pol, cor.plant = cor.plant))
	
			} else {
				return(c(NA,NA))
			}
		} else {
				return(c(NA,NA))
			}	
	
	} else {
		return(c(NA,NA))
	}	

}


##  bring everything together

master.fun <- function(tree1,tree2) {
	
	these.trait <- character.evo(tree1, tree2)
	
	same.same.stat <- link.rule(these.trait$same.same)
	same.diff.stat <- link.rule(these.trait$same.diff)
	diff.same.stat <- link.rule(these.trait$diff.same)
	diff.diff.stat <- link.rule(these.trait$diff.diff)

	return(cbind(same.same.stat,same.diff.stat,diff.same.stat,diff.diff.stat))
}




