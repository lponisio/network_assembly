
## interaction matrices are matrices of an arbitary dimension with
## cells indicating the number of interactions between the nth and mth
## species. These functions can handle binary, integer and decimal
## interactions.

calc.metric <- function(dat.web) {
    ## calculates greedy, random walk and edge betweenness modularly
    ## algoithums and NODF
    ## takes an interaction matrix, returns the four metrics (three
    ## modularity and one nestedness)

    ## calculates modularity
    calc.mod <- function(dat.web){
        ## calculates modularity

        mut.adj <- function(x) {
            ## converts a p-a matrix to a graph for modularity computation
            nr <- dim(x)[1]
            nc <- dim(x)[2]
            to.fill <- matrix(0, ncol=nc + nr, nrow=nc + nr)
            to.fill[1:nr,(nr+1):(nc+nr)] <- x
            adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
            return(adj.mat)
        }
        graph <- mut.adj(dat.web)
        weights <- as.vector(dat.web)
        weights <- weights[weights != 0]

        ## if matrix is binary, modularity calculate is not affected by
        ## the weights
        greedy <- modularity(graph,
                             membership(fastgreedy.community(graph,
                                                          weights=weights)),
                             weights=weights)

        random.walk <-  modularity(graph,
                                   membership(walktrap.community(graph,
                                                            weights=weights)),
                                   weights=weights)
        dendro <-  modularity(graph,
                              membership(edge.betweenness.community(graph,
                                                            weights=weights)),
                              weights=weights)
        return(c(greedy, random.walk, dendro))
    }

    dat.web <- as.matrix(empty(dat.web))
    ## matrix of all the same number, return 0s
    if(min(dat.web) == max(dat.web)){
        return(c(mets=rep(0, 1),
                 mod.met=rep(0,3)))
        ## if a binary (or only two values) matrix, use unweighted nodf
    }else if(min(dat.web[dat.web != 0]) ==
             max(dat.web[dat.web !=0])){
        mets <- nestednodf(dat.web,
                           weighted=FALSE)$statistic["NODF"]
        ## if a non binary matrix
    } else{
        mets <- nestednodf(dat.web,
                           weighted=TRUE)$statistic["NODF"]
    }
    mod.met <- calc.mod(dat.web)
    return(c(mets, mod.met= mod.met))

}

prob.null <- function(M) {
    ## null model for reshuffling comunities based on interaction
    ## frequencies. Modified from wine function in bipartite package
    randiag <- function(vdims) {
        mats <- diag(vdims)
        ## what cell number are the interactions
        posdiag <- 1:vdims + ((1:vdims - 1) * vdims)
        desp <- matrix(rep(1:vdims, vdims), vdims, vdims,
                       byrow = TRUE) - (1:vdims)
        fdesp <- sample(1:vdims)
        move <- desp[cbind(1:vdims, fdesp)]
        moved <- posdiag + move
        mdesp <- matrix(0, vdims, vdims)
        mdesp[moved] <- 1
        return(mdesp)
    }
    M <- as.matrix(M)
    values <- M[which(M > 0)]
    lvalues <- length(values)
    if (identical(dim(M)[1], dim(M)[2])) {
        vdims <- dim(M)[1]
        vdimb <- dim(M)[1]
    }
    if (!identical(dim(M)[1], dim(M)[2])) {
        dims <- which(dim(M) == min(dim(M)))
        vdims <- dim(M)[dims]
        dimb <- which(dim(M) == max(dim(M)))
        vdimb <- dim(M)[dimb]
    }
    MR <- matrix(0, vdims, vdimb)
    lMR <- vdims * vdimb
    sample1 <- sample(vdimb, vdims)
    diag1 <- randiag(vdims)
    MR[, sample1] <- diag1
    sample2 <- (1:vdimb)[-sample1]
    pos <- sample(vdims, length(sample2), replace = TRUE)
    MR[cbind(pos, sample2)] <- 2
    MRoccupied <- which(MR > 0)
    vleft <- lvalues - vdimb
    if (vleft > 0)
        MRoccupied <- c(MRoccupied, sample((1:lMR)[-MRoccupied],
                                           vleft))
    MR[MRoccupied] <- sample(values)
    if (dim(MR)[1] != dim(M)[1])
        MR <- t(MR)
    return(MR)
}


null.stat <- function(dat.web) {
    ## function to simulate 1 null, and calculate statistics on it
    ## takes interaction matrix, returns network metrics
    sim.web <- prob.null(dat.web)
    return(calc.metric(sim.web))
}

network.metrics <- function (dat.web, N, nmets=4) {
    ##  function that computes summary statistics on simulated null matrices
    ##  (nulls simulated from web N times)
    ## nmets is the number of metrics calcualted from the networks, in
    ##  this case it is always 4
    pvals <- function(stats, nnull){
        ## calculate pvalues
        rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
    }
    zvals <-function(stats){
        ## calculate zvalues based on SD
        z.sd <- (stats[,1] -
                 apply(stats, 1, mean, na.rm = TRUE))/
            apply(stats, 1, sd, na.rm = TRUE)
        z.sd[is.infinite(z.sd)] <- NA
        return(z.sd)
    }
    ## check that matrix is proper format (no empty row/col and no NAs)
    if(all(is.na(dat.web) == FALSE)) {
        ## drop empty rows and columns
        dat.web <- as.matrix(empty(dat.web))
        ## check to make sure emptied matrix is large enough
        ## to calculate statistics on
        if(is.matrix(dat.web)){
            if(all(dim(dat.web) >= 2)) {
                ## calculate null metrics
                null.stat <- replicate(N, null.stat(dat.web),
                                       simplify=TRUE)
                ## calculate metrics from data
                true.stat <- calc.metric(dat.web)
                out.mets <- cbind(true.stat, null.stat)
                ## compute z scores
                zvalues <- zvals(out.mets)
                ## compute p-values
                pvalues <- pvals(out.mets, N)
                return(cbind(true.stat, zvalues, pvalues))
            }
        }
    }
    ## 3 rows for the true metric, zvalue, and pvalue
    return(matrix(rep(NA,3*nmets), ncol=nmets))
}

