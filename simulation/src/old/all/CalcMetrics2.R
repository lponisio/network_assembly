## zscores are calculated including the true metric
## includes three methods for modularity computation, hierarchical
## clustering (edge.betweenness.community), dynamic algorithum
## (walktrap.community) and a greedy algorithum (fastgreedy.community)
## no longer calculate shannon div or evenness
## given inpute matrix (`data') returns desired summary statistics

calc.metric <- function(dat.web) {
  ## calculates modularity
  calc.mod <- function(dat.web){ 
    ## converts a p-a matrix to a graph for modularity computation 
    mut.adj <- function(x) {
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
                                                                weights=
                                                                weights)),
                          weights=weights)
    return(c(greedy, random.walk, dendro))
  }

  dat.web <- as.matrix(empty(dat.web))
  ## matrix of all the same number
  if(min(dat.web) == max(dat.web)){
    return(c(mets=rep(0, 1),
             mod.met=rep(0,3)))
    
  }else{
    ## wbinary=TRUE ensures binary matrices are calcualted like
    ## unweighted
    mets <- nestednodf(dat.web,
                       weighted=TRUE,
                       wbinary=TRUE)$statistic["NODF"]
  }
  mod.met <- calc.mod(dat.web)
  return(c(mets, mod.met= mod.met))
  
}

prob.null.binary <- function(probs, fill, nrow, ncol) {
  ones <- sample(1:length(probs), fill, replace=FALSE,
                 prob=probs)
  ## create resultant matrix
  interact <- matrix(0, nrow=nrow, ncol=ncol)
  interact[ones] <- 1
  return(interact)
}

prob.null.weighted <- function(probs, fill, num.int, nrow, ncol) {
  interact <- prob.null.binary(probs, num.int, nrow, ncol)
  int.probs <- probs[interact>0]
  int.probs <- int.probs/sum(int.probs)
  interact[interact > 0] <- int.probs*fill
  return(interact)
}

## function to simulate 1 null, and calculate statistics on it
null.stat <- function(prob, fill, num.int, nr, nc, binary) {
  if(binary){
    sim.web <- prob.null.binary(prob, fill, nr, nc)
    ## following is only necessary when using a null model that
    ## does not constraint marginals
    while(nrow(empty(sim.web)) < nr-5 |
          ncol(empty(sim.web)) < nc-5 |
          is.matrix(sim.web) == FALSE)
      sim.web <- prob.null.binary(prob, fill, nr, nc)
  } else {
    sim.web <- prob.null.weighted(prob, fill, num.int, nr, nc)
    while(nrow(empty(sim.web)) < nr-5 |
          ncol(empty(sim.web)) < nc-5 |
          is.matrix(sim.web) == FALSE)
      sim.web <- prob.null.weighted(prob, fill, nr, nc)
  }
  return(calc.metric(sim.web))
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
network.metrics <- function (dat.web, N) {
  ## calculate pvalues
  pvals <- function(stats, nnull){
    rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
  }
  ## calculate zvalues two different ways
  zvals <-function(stats){
    z.sd <- (stats[,1] -
             apply(stats, 1, mean, na.rm = TRUE))/
               apply(stats, 1, sd, na.rm = TRUE)
    z.sd[is.infinite(z.sd)] <- NA
    return(z.sd)
  }
  ## check that matrix is proper format (no empty row/col and no NAs)
  if(sum(dat.web) > 5 &
    all(is.na(dat.web) == FALSE)) {
    ## drop empty rows and columns
    dat.web <- as.matrix(empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(is.matrix(dat.web)){
      if(all(dim(dat.web) >= 5)) {
        ## produce N sets of statistics form N different nulls
        row.prob <- rowSums(dat.web)
        col.prob <- colSums(dat.web)
        dat.web.prob <- expand.grid(row.prob, col.prob)
        dat.web.prob <- apply(dat.web.prob, 1, mean)
        dat.web.fill <- sum(dat.web)
        dat.web.num.int <- sum(dat.web>0)
        ## calculate null metrics
        ifelse((min(dat.web[dat.web !=0]) == max(dat.web[dat.web !=0])),
               yes= binary <- TRUE,
               no= binary <- FALSE)
        null.stat <- replicate(N, null.stat(dat.web.prob,
                                            dat.web.fill,
                                            dat.web.num.int,
                                            nrow(dat.web),
                                            ncol(dat.web),
                                            binary),
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
  return(matrix(rep(NA,4*3), ncol=3)) 
}
