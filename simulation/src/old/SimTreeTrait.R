library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)

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
  tree2 <- reorder(tree2, "pruningwise")
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

interact.sample <- function(abundance, trait.overlap, np){
  trait <- matrix(trait.overlap, nrow=np)
  trait.abund <- trait*abundance
  return(list(evo = trait, eco = trait.abund))
}

##simulate abundances and define interaction probability 
sim.abund <- function(nanimal, nplant, mean.dist, sd.dist){
  abund.a <- rlnorm(nanimal, meanlog= mean.dist, sdlog= sd.dist)
  abund.p <- rlnorm(nplant, meanlog= mean.dist, sdlog= sd.dist)
  abund <- outer(abund.p, abund.a)
  return(abund)
}


matching.trait <- function(traits, combin){
  calc.range <- function(traits){
    ranges <-  cbind((traits - traits/2), (traits + traits/2))
    ranges[ranges < 0] <- 0
    return(ranges)
  }
  min.max <- function(prange, combin, index){
    pmin <- prange[combin[,index],1]
    pmax <- prange[combin[,index],2]
    return(rbind(pmin, pmax))
  } 

  ranges <- lapply(traits, FUN= calc.range)
  mms <- vector("list", length= 2)
  for (i in 1:2){
    mms[[i]] <- min.max(ranges[[i]], combin= combin, index= i)
  }

  lower.lim <- ifelse(mms[[1]][1,]  > mms[[2]][1,], mms[[1]][1,], mms[[2]][1,])
  upper.lim <- ifelse(mms[[1]][2,] < mms[[2]][2,], mms[[1]][2,], mms[[2]][2,])
  
  match.trait  <- upper.lim - lower.lim 
  match.trait[match.trait < 0] <- 0
  return(match.trait)
}

barrior.trait <- function(traits, combin){
  a.trait.combin <-traits[[1]][combin[,1]]
  p.trait.combin <- traits[[2]][combin[,2]]
  barrior <- ifelse(a.trait.combin > p.trait.combin, 1, 0)
  return(barrior)
}

##  takes two trait vectors (simulated from phylogenies) and computes
##  different plant-animal matrices based on different linkage rules

make.pa.mat <- function(traits, mean.abund, sd.abund){

  if(sum(unlist(traits)) != 0 &
     length(unlist(traits)[unlist(traits) != 0]) >= 5 &
     all(is.finite(unlist(traits)) == TRUE)){
    
    ns <- lapply(traits, FUN=length)

    abund.mat <- sim.abund(nanimal = ns[[1]], nplant = ns[[2]],
                           mean.dist = mean.abund, sd.dist = sd.abund)

    pa.combin <- expand.grid(1:ns[[1]] , 1:ns[[2]])

    evo.mats <- list(barrior = matrix(barrior.trait(traits,
                       pa.combin), ncol= ns[[2]]),
                     neutral = matrix(1, ncol = ns[[2]], nrow =
                       ns[[1]]),
                     matching = matrix(matching.trait(traits, pa.combin),
                       ncol = ns[[2]]))

    interact <-  lapply(evo.mats, interact.sample, abund = abund.mat, np = ns[[2]])

    return(unlist(interact, recursive = FALSE))

  }else{
    return(NA)
  }
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


##  uses `make.pa.mat' to produce a list of `fundamental' matrices
##  network statistics are then calculated on those matrices and returned

link.rule <- function(traits, nnull, meanAbund, sdAbund) {
  tree.dist <- traits[3:4]
  interact <- make.pa.mat(traits[1:2], meanAbund, sdAbund)
  ## calculate metrics for each matrix in fundamental
  out.metrics <- sapply(interact, network.metrics, N = nnull)

  rownames(out.metrics) <- c("NODF", "Shannon diversity", "interaction evenness", "H2",
                             "modularity", "z.NODF", "z.div", "z.even", "z.H2", "z.mod",
                             "z2.NODF", "z2.div", "z2.even", "z2.H2", "z2.mod", "p.NODF",
                             "p.div", "p.even", "p.H2", "p.mod")

  ## calculate phylo signal of interactions
  out.cors <- sapply(interact, function(x)
                     interact.phylo.signal(tree.dist,x))
  
  out.fill <- sapply(interact, sum)
  empty.interact <- lapply(interact, empty)
  out.npol <- sapply(empty.interact, ncol)
  out.nplant <- sapply(empty.interact, nrow)                   
  
  return(rbind(out.metrics, out.cors, out.fill, out.npol, out.nplant))
}

## bring everything together
master.fun <- function(tree1, tree2, nnull, mean.Abund, sd.Abund){
  these.trait <- character.evo(tree1, tree2)
  same.same.stat <- link.rule(these.trait$same.same, nnull, mean.Abund,
                              sd.Abund)
  same.diff.stat <- link.rule(these.trait$same.diff, nnull, mean.Abund,
                              sd.Abund)
  diff.same.stat <- link.rule(these.trait$diff.same, nnull, mean.Abund,
                              sd.Abund)
  diff.diff.stat <- link.rule(these.trait$diff.diff, nnull, mean.Abund,
                              sd.Abund)
  return(t(cbind(same.same.stat,
                 same.diff.stat,
                 diff.same.stat,
                 diff.diff.stat)))
}




