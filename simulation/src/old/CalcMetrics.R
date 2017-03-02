## added if statement to modularity calculation
## same as 'prob.null' but takes 'probs' and 'fill' as arguments to
## speed computation and returns quantitative matrix, instead of
## binary

library(bipartite)
library(igraph)

prob.null2 <- function(probs, fill, nrow, ncol) {
  interact <- rmultinom(1,fill,probs)[,1] 
  return(matrix(interact,nrow=nrow,ncol=ncol))
}

## converts a p-a matrix to a graph for modularity computation 
mut.adj <- function(x) {
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  
  to.fill <- matrix(0,ncol=nc + nr, nrow=nc + nr)
  
  to.fill[1:nr,(nr+1):(nc+nr)] <- x
  
  adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
  
  return(adj.mat)
}

## calculates modularity using edge-betweenness
calc.mod <- function(data){ 
  graph <- mut.adj(data)
  weights <- as.vector(data)
  weights <- weights[weights != 0]
  
  mod.met <- try(modularity(graph,
                            membership(edge.betweenness.community(graph)),
                            weights=weights), silent=TRUE)

  if(inherits(mod.met, "try-error")) {
    mod.met <- NA
  }
  return(mod.met)
}

## given inpute matrix (`data') returns desired summary statistics
calc.metric <- function(data) {
  data <- as.matrix(empty(data))
  
  network.metrics <- c("weighted NODF",
                       "Shannon diversity",
                       "interaction evenness", "H2")
  
  if(all(data == 1)){ # matrix of all ones return 0
    
    return(c(mets=rep(0, length(network.metrics)),
             mod.met=0))
    
  }else if(all(data[data !=0] == 1)){ # matrix all integers
    
    nodf <- nestednodf(data)$statistic["NODF"]
    mets <- c(nodf, try(networklevel(data,
                                     index=network.metrics[-1], H2_integer=TRUE), silent=TRUE))
    
  }else{
    
    mets <- try(networklevel(data, 
                             index=network.metrics, H2_integer=FALSE), silent=TRUE)
  }        

  if(inherits(mets, "try-error")) {
    mets <- rep(NA, length(network.metrics))
  }

  mod.met <- calc.mod(data)

  return(c(mets, mod.met=mod.met))
}


## function to simulate 1 null, and calculate statistics on it
null.stat <- function(prob, fill, nr, nc) {

  
  sim.web <- prob.null2(prob,fill,nr,nc)  # simulate null web, previous versions used r2dtable or vaznull.fast

                                        # following line is only necessary when using a null model that
                                        # does not constraint marginals, i.e. prob.null2
  while(any(dim(sim.web) < 5) | is.matrix(sim.web) == FALSE)
    sim.web <- prob.null2(prob,fill,nr,nc) 

  return(calc.metric(sim.web))
}

## calculate pvalues
pvals <- function(stats, nnull){
  rowSums(stats[,-1] >= stats[,rep(1, ncol(stats) -1)])/(nnull + 1)
}

## calculate zvalues two different ways
zvals <-function(stats){
  z.sd <- (stats[,1] -
           mean(stats[,-1],na.rm = TRUE))/sd(stats[,-1],na.rm = TRUE)
  
  z.mean <- (stats[,1] -
             mean(stats[,-1],na.rm = TRUE))/mean(stats[,-1],na.rm = TRUE)
  
  return(cbind(z.sd, z.mean))
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)

network.metrics <- function (data, N) {
                                        # check that matrix is proper format (no empty row/col and no NAs)
  if(sum(data != 0) & all(is.na(data) == FALSE)) {

                                        # drop empty rows and columns
    data <- as.matrix(empty(data))
                                        # check to make sure the rounded,emptied matrix is large enough
                                        # to calculate statistics on
    if(is.matrix(data)){
      if(all(dim(data) >= 5) ) {
                                        #  produce N sets of statistics form N different nulls
                                        #  should return matrix that is 5 row by N col
        row.prob <- rowSums(data)
        col.prob <- colSums(data)
        data.prob <- expand.grid(row.prob,col.prob)
        data.prob <- data.prob[,1]*data.prob[,2]
        data.fill <- sum(data)
        null.stat <- replicate(N, null.stat(data.prob, data.fill, # calculate null metrics
                                            nrow(data),ncol(data)),
                               simplify=TRUE)  
        
        true.stat <- calc.metric(data) # calculate metrics from data

        out.mets <- cbind(true.stat, null.stat)
        
        zvalues <- zvals(out.mets) # compute z scores
        
        pvalues <- pvals(out.mets, N)  # compute p-values

        return(cbind(true.stat, zvalues, pvalues))
      }
    }
  }
  blank <- matrix(rep(NA,5*4), ncol=4)
  colnames(blank) <- c("true.stat", "z.sd", "z.mean", "pvalues")
  rownames(blank) <-  c("weighted NODF",
                       "Shannon diversity",
                       "interaction evenness", "H2", "mod.met")
  return(blank) 
}
