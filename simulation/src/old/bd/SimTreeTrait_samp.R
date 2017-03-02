## simulate trees
sim.phylo <- function(nsim, mu, lambda, age, nspecies) {
  
  tree.1 <- sim.bd.taxa.age(n=nspecies, numbsim=nsim, lambda=lambda,
                            mu= mu, frac= 1, age= age, mrca= FALSE)
  
  tree.2 <- sim.bd.taxa.age(n=nspecies, numbsim=nsim, lambda=lambda,
                            mu= mu, frac= 1, age= age, mrca= FALSE) 

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
                                          model= "REML")$ace), silent= TRUE)
  
  while(inherits(optima.tree2, "try-error")) {
    print("crazy tree!")
    sample.trait <- sample(tree1.trait1, size= Ntip(tree2), replace= TRUE)
    names(sample.trait) <- tree2$tip.label
    tree2 <- reorder(tree2, "pruningwise")
    optima.tree2 <- try(c(sample.trait, ace(x= sample.trait, phy=
                                            tree2, model=
                                            "REML")$ace), silent= TRUE)
  }

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

make.pa.mat <- function(prms, traits, mean.abund, sd.abund, w.evo){
  ## takes overlap matrices and abundance matrix and calculates the
  ## overlap*abundance and sampled matrices
  interact.sample <- function(abundance, trait.overlap, nplant, w.evo){
    trait <- matrix(trait.overlap, nrow=nplant)
    trait <- (trait/sum(trait))*(100*w.evo) 
    abundance <- (abundance/sum(abundance))*100
    ## ecological matrix with random abundances
    trait.abund <- trait*abundance
    trait.abund <-  (trait.abund/sum(trait.abund))*100

    ## ecological matrix with abundance correlated with traits
    trait.2 <- trait[order(rowSums(trait), decreasing=TRUE),]
    trait.2 <- trait.2[,order(colSums(trait.2), decreasing =TRUE)]

    abund.2 <- abundance[order(rowSums(abundance), decreasing=TRUE),]
    abund.2 <- abund.2[,order(colSums(abund.2), decreasing =TRUE)]

    trait.abund2 <- trait.2*abund.2
    trait.abund2 <-  (trait.abund2/sum(trait.abund2))*100

    ##sampling

    fill.trait.abund  <- sum(trait.abund, na.rm=TRUE)
    fill.trait.abund2 <- sum(trait.abund2, na.rm=TRUE)

       if(fill.trait.abund !=0 & is.finite(fill.trait.abund)){
         trait.abund.samp.2 <- matrix(rmultinom(1,
                                                fill.trait.abund*0.2,
                                                trait.abund),
                                      nrow=nplant)
         trait.abund.samp.4 <- matrix(rmultinom(1,
                                                fill.trait.abund*0.4,
                                                trait.abund),
                                      nrow=nplant)
         trait.abund.samp.6 <- matrix(rmultinom(1,
                                                fill.trait.abund*0.6,
                                                trait.abund),
                                      nrow=nplant)
         trait.abund.samp.8 <- matrix(rmultinom(1,
                                                fill.trait.abund*0.8,
                                                trait.abund),
                                      nrow=nplant)
         
         trait.abund2.samp.2 <- matrix(rmultinom(1,
                                                 fill.trait.abund2*0.2,
                                                 trait.abund2),
                                       nrow=nplant)
         trait.abund2.samp.4 <- matrix(rmultinom(1,
                                                 fill.trait.abund2*0.4,
                                                 trait.abund2),
                                       nrow=nplant)
         trait.abund2.samp.6 <- matrix(rmultinom(1,
                                                 fill.trait.abund2*0.6,
                                                 trait.abund2),
                                       nrow=nplant)
         trait.abund2.samp.8 <- matrix(rmultinom(1,
                                                 fill.trait.abund2*0.8,
                                                 trait.abund2),
                                       nrow=nplant)
         
         return(list(abund.samp.2 = trait.abund.samp.2,
                     abund.samp.4 = trait.abund.samp.4,
                     abund.samp.6 = trait.abund.samp.6,
                     abund.samp.8 = trait.abund.samp.8,
                     abund2.samp.2 = trait.abund2.samp.2,
                     abund2.samp.4 = trait.abund2.samp.4,
                     abund2.samp.6 = trait.abund2.samp.6,
                     abund2.samp.8 = trait.abund2.samp.8))
       } else{
         return(list(abund.samp.2 = 0,
                     abund.samp.4 = 0,
                     abund.samp.6 = 0,
                     abund.samp.8 = 0,
                     abund2.samp.2 = 0,
                     abund2.samp.4 = 0,
                     abund2.samp.6 = 0,
                     abund2.samp.8 = 0))
       }
     }

    ##simulate abundances and define interaction probability 
    sim.abund <- function(nanimal, nplant, mean.dist, sd.dist){
      abund.a <- rpoilog(S = nanimal, mu = mean.dist, sig= sd.dist,
                         keep0 =TRUE)
      abund.p <- rpoilog(S = nplant, mu = mean.dist, sig= sd.dist,
                         keep0=TRUE)
      abund <- outer(abund.p, abund.a) ## plants are the rows
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
       all(is.finite(unlist(traits)))){
      
      abund.mat <- sim.abund(nanimal = prms$sp, nplant = prms$sp, 
                             mean.dist = mean.abund, sd.dist = sd.abund)

      evo.mats <- list(barrior= barrior.trait(prms, traits),
                       neutral= matrix(1, ncol=prms$sp, nrow=prms$sp),
                       matching= matching.trait(prms, traits))
      
      interact <- lapply(evo.mats, interact.sample, abund = abund.mat,
                         nplant = prms$sp, w.evo= w.evo)
      
      return(unlist(interact, recursive = FALSE))

    }
    return(NA)
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

  link.rule <- function(prms, traits, nnull, meanAbund, sdAbund, w.evo) {
    tree.dist <- traits[3:4]
    interact <- make.pa.mat(prms, traits[1:2], meanAbund, sdAbund, w.evo)
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
  master.fun <- function(prms, tree1, tree2, nnull, mean.Abund,
                         sd.Abund, w.evo) {
    these.trait <- character.evo(tree1, tree2)
    same.same.stat <- link.rule(prms,
                                these.trait$same.same, nnull,
                                mean.Abund, sd.Abund, w.evo)
    same.diff.stat <- link.rule(prms,
                                these.trait$same.diff, nnull,
                                mean.Abund, sd.Abund, w.evo)
    diff.same.stat <- link.rule(prms,
                                these.trait$diff.same, nnull,
                                mean.Abund, sd.Abund, w.evo)
    diff.diff.stat <- link.rule(prms,
                                these.trait$diff.diff, nnull,
                                mean.Abund, sd.Abund, w.evo)
    
    return(t(cbind(same.same.stat,
                   same.diff.stat,
                   diff.same.stat,
                   diff.diff.stat))) 
  }
