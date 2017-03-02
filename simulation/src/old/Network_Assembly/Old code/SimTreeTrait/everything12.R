## avoid calculating metrics on neutral fund matrix
## let trees evolve and then drop tips
# nnull added
##changed the way different-same traits are simulated
## no sampling

## needed libraries
library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)
library(distr)


## simulate trees
sim.phylo <- function(nsim, mu, lambda, age, nspecies){
	
	tree.1 <- sim.bd.taxa.age(n=1000, numbsim=nsim, lambda=lambda,
                            mu=mu, frac=1, age=age, mrca=TRUE)                     
  
  tree.2 <- sim.bd.taxa.age(n=1000, numbsim=nsim, lambda=lambda,
                            mu=mu, frac=1, age=age, mrca=TRUE) 
                            
    tree.1.trim <- lapply(tree.1, FUN=function(X) drop.tip(X, sample(X$tip.label, 1000-nspecies)))
     
    tree.2.trim <- lapply(tree.2, FUN=function(X) drop.tip(X, sample(X$tip.label, 1000-nspecies)))
     
  return(list(tree.1.trim, tree.2.trim))
}

## simulate character evolution
character.evo <- function(tree1, tree2){
  tree1.trait1 <- rTraitCont(tree1, model="BM")
  tree1.trait2 <- rTraitCont(tree1, model="BM")
  tree2.trait1 <- rTraitCont(tree2, model="BM")
 
  sample.trait <- sample(tree1.trait1, size=Ntip(tree2), replace=TRUE)
  names(sample.trait) <- tree2$tip.label
  tree2 <- reorder(tree2, "p")
  optima.tree2 <- c(sample.trait,ace(x=sample.trait, phy=tree2, model="REML")$ace)
  optima.tree2 <- optima.tree2[tree2$edge[,2]]
  tree2.trait2 <- rTraitCont(tree2, model="OU", theta=optima.tree2)
  
  ## phylogenetic distances between tips in plant and pol phylos
  tree1.dist <- cophenetic(tree1)
  tree2.dist <- cophenetic(tree2)
  
  same.same <- list(tree1.trait1, tree1.trait1, tree1.dist, tree1.dist)
  same.diff <- list(tree1.trait1, tree1.trait2, tree1.dist, tree1.dist)		
  diff.diff <- list(tree1.trait1, tree2.trait1, tree1.dist, tree2.dist)
  diff.same <- list(tree1.trait1, tree2.trait2, tree1.dist, tree2.dist)
  
  out.traits <- list(same.same=same.same,
              	same.diff=same.diff, 
              	diff.diff=diff.diff,
              	diff.same=diff.same)
    	
  return(out.traits)
}

## takes overlap matrices and abundance matrix and calculates the
## overlap*abundance and sampled matrices

interact.sample <- function(abundance, trait.overlap, npol, nplant){
  
  trait <- matrix(trait.overlap, nrow=npol)
  
  if(sum(trait) != 0 & length(trait[trait != 0]) >= 5 & all(is.finite(trait) == TRUE)){  
    
    trait.abund <- trait*abundance
 	return(list(trait, trait.abund))
    
  } else{
    return(list(matrix(0, nrow=npol, ncol=nplant), matrix(0, nrow=npol, ncol=nplant)))
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
  
                                        # narrow
  range.pol.narrow <- cbind((pol.trait - pol.trait/2), (pol.trait + pol.trait/2))
  range.pol.narrow[range.pol.narrow < 0] <- 10^-10
  
  range.plant.narrow <-  cbind((plant.trait - plant.trait/2), (plant.trait + plant.trait/2))
  range.plant.narrow[range.plant.narrow < 0] <- 10^-10
  
  ## narrow traits 
  pol.narrow.min <- range.pol.narrow[pa.comb[,1],1]
  pol.narrow.max <- range.pol.narrow[pa.comb[,1],2]
  
  plant.narrow.min <- range.plant.narrow[pa.comb[,2],1]
  plant.narrow.max <- range.plant.narrow[pa.comb[,2],2]

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
  
  ## continuous with trait barrior
  
  plant.trait.combin <- plant.trait[pa.comb[,2]]
  pol.trait.combin <- plant.trait[pa.comb[,1]]
  
  barrior <- ifelse(pol.trait.combin > plant.trait.combin, 1, 0)
  
  trait.bar <- interact.sample(abund, barrior, npol, nplant)
  
  neutral <- interact.sample(abund, matrix(1, ncol=nplant, nrow=npol), npol, nplant)

  interact <- c( 
	 	trait.bar=trait.bar,
	 	neutral=neutral, 
	 	trait.narrow=trait.narrow
	 	 )
  return(interact)
}

##  uses `make.pa.mat' to produce a list of `fundamental' matrices
##  network statistics are then calculated on those matrices and returned

link.rule <- function(traits, nnull) {
  tree.dist <- traits[3:4]
  interact <- make.pa.mat(traits[1:2])
  
  ## calculate metrics for each matrix in fundamental
  out.metrics <- sapply(interact, network.metrics, N=prms$nnull)

  ## calculate phylo signal of interactions
  out.cors <- sapply(interact, function(x)
                     interact.phylo.signal(tree.dist,x))
                     
  out.fill <- sapply(interact, sum)
	
	out.npol <- sapply(interact, ncol)
	
	out.nplant <- sapply(interact, nrow)                   

  return(rbind(out.metrics, out.cors, out.fill, out.npol, out.nplant))
}

## calculate phylogenetic structure of interactions tree 1=polinators,
## tree 2=plants, int.mat is a list of interaction matrices

interact.phylo.signal <- function(dist.tree, int.mat) {
  dist.tree1 <- dist.tree[[1]]
  dist.tree2 <- dist.tree[[2]]
  
  if(sum(int.mat) != 0 & all(is.na(int.mat) == FALSE)){
    
    cNonEmpty <- colSums(int.mat) > 0
    rNonEmpty <- rowSums(int.mat) > 0
    
    int.mat <- int.mat[rNonEmpty,cNonEmpty]
    dist.tree1 <- dist.tree1[cNonEmpty,cNonEmpty]
    dist.tree2 <- dist.tree2[rNonEmpty, rNonEmpty]
    
    if(is.matrix(int.mat)){
      
      if(all(dim(int.mat) >= 10) & any(int.mat!=1)) {

        plant.dist <- as.matrix(vegdist(int.mat, method="bray",
                                        diag=TRUE, upper=TRUE),
                                		ncol=ncol(dist.tree2))
        pol.dist <- as.matrix(vegdist(t(int.mat), method="bray",
                                      diag=TRUE, upper=TRUE),
                                		ncol=ncol(dist.tree1))
	options(warn =-1)
        cor.plant <- mantel(plant.dist, dist.tree2,
                            permutations=FALSE,
                            method="spearman")$statistic
        cor.pol <- mantel(pol.dist, dist.tree1,
                          permutations=FALSE,
                          method="spearman")$statistic 
	
        return(c(cor.pol=cor.pol, cor.plant=cor.plant))
      }
    }
  }
  return(c(NA,NA))
}

## bring everything together
master.fun <- function(tree1,tree2) {
  these.trait <- character.evo(tree1, tree2)
  same.same.stat <- link.rule(these.trait$same.same)
  same.diff.stat <- link.rule(these.trait$same.diff)
  diff.same.stat <- link.rule(these.trait$diff.same)
  diff.diff.stat <- link.rule(these.trait$diff.diff)
  return(t(cbind(same.same.stat,
                 same.diff.stat,
                 diff.same.stat,
                 diff.diff.stat)))
}




