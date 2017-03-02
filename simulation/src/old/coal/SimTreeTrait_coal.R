## simulate trees
sim.phylo <- function(nspecies) {
  tree.1 <- rcoal(n= nspecies)
  tree.2 <- rcoal(n= nspecies)
  return(list(tree.1, tree.2))
}

## simulate character evolution
character.evo <- function(tree1, tree2){
  tree1.trait1 <- rTraitCont(tree1, model="BM")
  tree1.trait2 <- rTraitCont(tree1, model="BM")
  tree2.trait1 <- rTraitCont(tree2, model="BM")
  
  sample.trait <- sample(tree1.trait1, size= Ntip(tree2), replace= TRUE)
  names(sample.trait) <- tree2$tip.label
  tree2 <- reorder(tree2, "pruningwise")
  
  optima.tree2 <- try(c(sample.trait, ace(x= sample.trait, phy= tree2,
                                          model= "REML")$ace), silent=
                                          TRUE)
  
  while(inherits(optima.tree2, "try-error")) {
    print("crazy tree!")
    sample.trait <- sample(tree1.trait1, size= Ntip(tree2), replace= TRUE)
    names(sample.trait) <- tree2$tip.label
    tree2 <- reorder(tree2, "pruningwise")
    optima.tree2 <- try(c(sample.trait, ace(x= sample.trait, phy=
                                            tree2, model=
                                            "REML")$ace), silent=
                                            TRUE) }
  
  optima.tree2 <- optima.tree2[tree2$edge[,2]]
  tree2.trait2 <- rTraitCont(tree2, model="OU", theta= optima.tree2)
  
  ## phylogenetic distances between tips in plant and pol phylos
  tree1.dist <- cophenetic(tree1)
  tree2.dist <- cophenetic(tree2)
  
  same.same <- list(tree1.trait1, tree1.trait1, tree1.dist, tree1.dist)
  same.diff <- list(tree1.trait1, tree1.trait2, tree1.dist, tree1.dist)		
  diff.diff <- list(tree1.trait1, tree2.trait1, tree1.dist, tree2.dist)
  diff.same <- list(tree1.trait1, tree2.trait2, tree1.dist, tree2.dist)
  
  out.traits <- list(same.same= same.same,
                     same.diff= same.diff, 
                     diff.diff= diff.diff,
                     diff.same= diff.same)

  return(out.traits)
}


##  takes two trait vectors (simulated from phylogenies) and computes
##  different plant-animal matrices based on different linkage rules

make.pa.mat <- function(prms, traits, mean.abund, sd.abund){

  ## takes overlap matrices and abundance matrix and calculates the
  ## overlap*abundance and sampled matrices
  interact.sample <- function(abundance, trait.overlap, np){
    trait <- matrix(trait.overlap, nrow=np)
    trait.abund <- trait*abundance
    
    return(list(evo = trait, eco = trait.abund))
  }

  ##simulate abundances and define interaction probability 
  sim.abund <- function(nanimal, nplant, mean.dist, sd.dist){
    abund.a <- rlnorm(nanimal, meanlog=mean.dist, sdlog= sd.dist)
    abund.p <- rlnorm(nplant, meanlog=mean.dist, sdlog= sd.dist)
    abund <- outer(abund.p, abund.a)
    return(abund)
  }

  matching.trait <- function(prms, traits){
    calc.range <- function(traits)
      cbind((traits - traits/2), (traits + traits/2))
    ranges <- lapply(traits, calc.range)
    ## the above assumption needs to be thought about
    plants  <- ranges[[1]][prms$combinations[,1],]
    animals <- ranges[[2]][prms$combinations[,2],]

    ## computes upper or lower limits, depending on f1 and f2
    lims <- function(f1, f2)
      mapply(function(a,b) f1(a,b),
             apply(plants, 1, f2), apply(animals, 1, f2))

    match.trait <- lims(max, min) - lims(min, max)
    match.trait[match.trait < 0] <- 0
    return(matrix(match.trait, prms$sp, prms$sp))
  }

  barrior.trait <- function(prms, traits){
    a.trait.combin <- traits[[1]][prms$combinations[,1]]
    p.trait.combin <- traits[[2]][prms$combinations[,2]]
    return(matrix(a.trait.combin > p.trait.combin,
                  nrow=prms$sp, ncol=prms$sp)*1)
  }
  
  if(sum(unlist(traits)) != 0 &
     length(unlist(traits)[unlist(traits) != 0]) >= 5 &
     all(is.finite(unlist(traits)) == TRUE)){
    
    abund.mat <- sim.abund(nanimal = prms$sp, nplant = prms$sp, 
                           mean.dist = mean.abund, sd.dist = sd.abund)

    evo.mats <- list(barrior= barrior.trait(prms, traits),
                     neutral= matrix(1, ncol=prms$sp, nrow=prms$sp),
                     matching= matching.trait(prms, traits))

    interact <-  lapply(evo.mats, interact.sample, abund = abund.mat,
                        np = prms$sp) 

    return(unlist(interact, recursive = FALSE))

  } else{
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

link.rule <- function(prms, traits, nnull, meanAbund, sdAbund) {
  tree.dist <- traits[3:4]
  interact <- make.pa.mat(prms, traits[1:2], meanAbund, sdAbund)
  ## calculate metrics for each matrix in fundamental
  out.metrics <- sapply(interact, network.metrics, N = nnull)

  rownames(out.metrics) <- c("NODF", "Shannon diversity",
                             "interaction evenness", "H2",
                             "modularity", "z.NODF", "z.div",
                             "z.even", "z.H2", "z.mod", "z2.NODF",
                             "z2.div", "z2.even", "z2.H2", "z2.mod",
                             "p.NODF", "p.div", "p.even",
                             "p.H2", "p.mod")

  ## calculate phylo signal of interactions
  out.cors <- sapply(interact, function(x)
                     interact.phylo.signal(tree.dist,x))
  
  out.fill <- sapply(interact, sum)
  empty.interact <- lapply(interact, empty)
  out.npol <- sapply(empty.interact, ncol)
  out.nplant <- sapply(empty.interact, nrow)                   

  return(rbind(out.metrics, out.cors, out.fill,
               out.npol, out.nplant))
}

## bring everything together
master.fun <- function(prms, tree1, tree2, nnull, mean.Abund, sd.Abund) {
  these.trait <- character.evo(tree1, tree2)
  same.same.stat <- link.rule(prms,
                              these.trait$same.same, nnull,
                              mean.Abund, sd.Abund)
  same.diff.stat <- link.rule(prms,
                              these.trait$same.diff, nnull,
                              mean.Abund, sd.Abund)
  diff.same.stat <- link.rule(prms,
                              these.trait$diff.same, nnull,
                              mean.Abund, sd.Abund)
  diff.diff.stat <- link.rule(prms,
                              these.trait$diff.diff, nnull,
                              mean.Abund, sd.Abund)
                              
  return(t(cbind(same.same.stat,
                 same.diff.stat,
                 diff.same.stat,
                 diff.diff.stat))) 
}
