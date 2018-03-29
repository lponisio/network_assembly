
character.evo <- function(tree1, tree2, sigma, alpha=1){
    ## simulate character evolution
    ## function takes two trees, and the variance of brownian motion
    tree1.trait1 <- rTraitCont(tree1, model='BM', sigma=sigma)
    tree1.trait2 <- rTraitCont(tree1, model='BM', sigma=sigma)
    tree2.trait1 <- rTraitCont(tree2, model='BM', sigma=sigma)

    ou.traits <- function(tree1.trait1, tree2){
        ## function to generate traits based on an Ornstein Uhlenbeck
        ## process
        gen.optima <- function(){
            sample.trait <- sample(tree1.trait1, size=Ntip(tree2),
                                   replace=TRUE)
            names(sample.trait) <- tree2$tip.label
            tree2 <- reorder(tree2, 'pruningwise')
            optima.tree2 <- c(sample.trait, ace(x=sample.trait,
                                                phy=tree2,
                                                method='pic')$ace)
            return(optima.tree2)
        }
        optima.tree2 <- gen.optima()
        optima.tree2 <- optima.tree2[tree2$edge[,2]]
        tree2.trait2 <- rTraitCont(tree2, model='OU',
                                   theta=optima.tree2, alpha=alpha)
        return(tree2.trait2)
    }
    tree2.trait2 <- ou.traits(tree1.trait1, tree2)
    tree1.trait1.2 <-  ou.traits(tree2.trait2, tree1)
    ## phylogenetic distances between tips in plant and pol phylos
    tree1.dist <- cophenetic(tree1)
    tree2.dist <- cophenetic(tree2)

    same.same <- list(tree1.trait1,
                      tree1.trait1,
                      tree1.dist,
                      tree1.dist)
    same.diff <- list(tree1.trait1,
                      tree1.trait2,
                      tree1.dist,
                      tree1.dist)
    diff.diff <- list(tree1.trait1,
                      tree2.trait1,
                      tree1.dist,
                      tree2.dist)
    diff.same <- list(tree1.trait1.2,
                      tree2.trait2,
                      tree1.dist,
                      tree2.dist)

    out.traits <- list(same.same=same.same,
                       same.diff=same.diff,
                       diff.diff=diff.diff,
                       diff.same=diff.same)
    return(out.traits)
}


make.pa.mat <- function(prms, traits){
    ## takes two trait vectors (simulated from phylogenies) and
    ## computes different plant-animal matrices based on different
    ## linkage rules

    calc.range <- function(trait) {
        cbind((trait - abs(trait*prms$range.size)),
        (trait + abs(trait*prms$range.size)))
    }
    ranges <- lapply(traits, calc.range)
    plants  <- ranges[[1]][prms$combinations[,1],]
    animals <- ranges[[2]][prms$combinations[,2],]

    ## computes upper or lower limits, depending on f1 and f2
    ## range is the min of the maxes - the max of the mins
    lims <- function(f1, f2)
        mapply(f1, apply(plants, 1, f2), apply(animals, 1, f2))

    match.trait <- lims(min, max) - lims(max, min)
    ## traits do not overlap if the above is negative
    match.trait[match.trait < 0] <- 0
    match.trait <- round(match.trait*1000, 0)
    ## only binary
    match.trait.bin <- match.trait
    match.trait.bin[match.trait.bin > 0] <- 1
    return(list(quan=matrix(match.trait, prms$sp, prms$sp),
                qual=matrix(match.trait.bin, prms$sp, prms$sp)))
}


interact.phylo.signal <- function(dist.tree, int.mat,
                                  dis.method="gower",
                                  cor.method="spearman") {
    ## calculate phylogenetic structure of interactions tree
    ## 1=polinators, tree 2=plants, int.mat is a list of interaction
    ## matrices

    dist.tree1 <- dist.tree[[1]]
    dist.tree2 <- dist.tree[[2]]

    if(sum(int.mat)){

        cNonEmpty <- colSums(int.mat) > 0
        rNonEmpty <- rowSums(int.mat) > 0

        int.mat <- empty(int.mat)
        dist.tree1 <- dist.tree1[cNonEmpty, cNonEmpty]
        dist.tree2 <- dist.tree2[rNonEmpty, rNonEmpty]

        if(is.matrix(int.mat)){
            if(all(dim(int.mat) >=5)) {
                plant.dist <- as.matrix(vegdist(int.mat,
                                                method=dis.method,
                                                diag=TRUE,
                                                upper=TRUE),
                                        ncol=ncol(dist.tree2))
                pol.dist <- as.matrix(vegdist(t(int.mat),
                                              method=dis.method,
                                              diag=TRUE,
                                              upper=TRUE),
                                      ncol=ncol(dist.tree1))
                options(warn=-1)
                cor.plant <- mantel(plant.dist,
                                    dist.tree2,
                                    permutations=0,
                                    method=cor.method)$statistic
                cor.pol <- mantel(pol.dist,
                                  dist.tree1,
                                  permutations=0,
                                  method=cor.method)$statistic
                ## subtrack 1 from dissimilarity to get simmilarity
                return(c(cor.pol=cor.pol,
                         cor.plant=cor.plant,
                         mean.dis.pol=1-mean(pol.dist),
                         mean.dis.plant=1- mean(plant.dist)))
            }
        }
    }
    return(c(NA, NA, NA, NA))
}


link.rule <- function(prms,
                      traits,
                      nnull) {
    ##  uses `make.pa.mat' to produce a list of `fundamental' matrices
    ##  network statistics are then calculated on those matrices and
    ##  returned
    tree.dist <- traits[3:4]
    interact <- make.pa.mat(prms,
                            traits[1:2])
    ## calculate metrics for each matrix in fundamental
    out.metrics <- sapply(interact,
                          network.metrics,
                          N=nnull)

    ## calculate phylo signal of interactions
    out.cors <- sapply(interact, function(x)
        interact.phylo.signal(tree.dist,x))
    out.fill <- sapply(interact, sum)
    empty.interact <- lapply(interact, empty)
    ## number of higher level species
    out.npol <- sapply(empty.interact, ncol)
    ## number of lower level species
    out.nplant <- sapply(empty.interact, nrow)
    ## mean degree higher level
    mean.int.pol <- rep(mean(apply(empty.interact$qual, 2, mean)), 2)
    ## mean degree lower level
    mean.int.plant <- rep(mean(apply(empty.interact$qual, 1, mean)),
                          2)
    ## connectance
    out.connect <- sapply(empty.interact, function(x){
        length(x[x !=0])/(nrow(x)*ncol(x))
    })
    return(rbind(out.metrics,
                 out.cors,
                 out.fill,
                 out.npol,
                 out.nplant,
                 mean.int.pol,
                 mean.int.plant,
                 out.connect))
}

master.fun <- function(prms,
                       tree1,
                       tree2,
                       nnull) {
    ## bring everything together
    these.trait <- character.evo(tree1, tree2, prms$sigma)
    same.same.stat <- link.rule(prms,
                                these.trait$same.same,
                                nnull)
    same.diff.stat <- link.rule(prms,
                                these.trait$same.diff,
                                nnull)
    diff.same.stat <- link.rule(prms,
                                these.trait$diff.same,
                                nnull)
    diff.diff.stat <- link.rule(prms,
                                these.trait$diff.diff,
                                nnull)
    return(t(cbind(same.same=same.same.stat,
                   same.diff=same.diff.stat,
                   diff.same=diff.same.stat,
                   diff.diff=diff.diff.stat)))
}
