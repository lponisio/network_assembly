## same as 'prob.null' but takes 'probs' and 'fill' as arguments to
## speed computation and returns quantitative matrix, instead of
## binary

prob.null2 <- function(probs,fill,nrow,ncol) {
  
  ##  draw random interactions
  interact <- rmultinom(1,fill,probs)[,1]
  return(matrix(interact,nrow=nrow,ncol=ncol))
}

## converts a p-a matrix to a graph for modularity computation 

mut.adj <- function(x) { ##and p-a matrix
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  
  to.fill <- matrix(0,ncol=nc + nr, nrow=nc + nr)
  
  to.fill[1:nr,(nr+1):(nc+nr)] <- x
  
  adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
  
  return(adj.mat)
}

## given inpute matrix (`data') returns desired summary statistics

calc.metric <- function(data) {
  
  data <- as.matrix(empty(data))
  
  mets <- networklevel(data, index=c("weighted NODF", "H2", "ISA"))
  
  graph <- mut.adj(data)
  weights <- as.vector(data)
  weights <- weights[weights != 0]
  
  eb <- edge.betweenness.community(graph)
  memb <- community.to.membership(graph, eb$merges, steps= (nrow(eb$merges) -1) )
  mod.met <- modularity(graph, memb$membership, weights=weights)
  
  return(c(mets, mod.met=mod.met))
}


## function to simulate 1 null, and calculate statistics on it
null.stat <- function(prob,fill,nr,nc) {
  
  ## simulate null web, previous versions used r2dtable or vaznull.fast
  sim.web <- prob.null2(prob,fill,nr,nc)
   
  ## following line is only necessary when using a null model that
  ## does not constraint marginals, i.e. prob.null2
  while(any(dim(sim.web) < 5) | is.matrix(sim.web) == FALSE)
    sim.web <- prob.null2(prob,fill,nr,nc) 
  
  return(calc.metric(sim.web))
}


##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)

network.metrics <- function (data, N) {
  
  ## print(data)

  ## check that matrix is proper format (no empty row/col and no NAs)
  if(sum(data != 0) & all(is.na(data) == FALSE)) {

    ## drop empty rows and columns
    data <- as.matrix(empty(data))
    
    ## check to make sure the rounded,emptied matrix is large enough
    ## to calculate statistics on 
    
    if(is.matrix(data)){
      if(all(dim(data) >= 5) ) {
        
        ## print(data)
    	
        ##  produce N sets of statistics form N different nulls
        ##  should return matrix that is 5 row by N col
        
        row.prob <- rowSums(data)
        col.prob <- colSums(data)
        data.prob <- expand.grid(row.prob,col.prob)
        data.prob <- data.prob[,1]*data.prob[,2]
        
        data.fill <- sum(data)
        
        null.stat <- replicate(N, null.stat(data.prob, data.fill,
                                            nrow(data),ncol(data)),
                               simplify=TRUE)  
##         print("null_met")
        true.stat <- calc.metric(data)
##         print("true_met")
        ##  combind into one matrix from which p-values will be calculated
        out.mets <- cbind(true.stat,null.stat)
        
        ##  compute z scores
        z.nodf <- (out.mets[1,1] - mean(out.mets[1,-1],na.rm=TRUE))/sd(out.mets[1,-1],na.rm=TRUE)
        z.h2 <- (out.mets[2,1] - mean(out.mets[2,-1],na.rm=TRUE))/sd(out.mets[2,-1],na.rm=TRUE)
        z.ISA <- (out.mets[3,1] - mean(out.mets[3,-1],na.rm=TRUE))/sd(out.mets[3,-1],na.rm=TRUE)       
        z.mod <- (out.mets[4,1] - mean(out.mets[4,-1],na.rm=TRUE))/sd(out.mets[4,-1],na.rm=TRUE)
        
        ##compute p-values
        pval.nodf <- sum(out.mets[1,-1] >= out.mets[1,1], na.rm=TRUE)/(N+1)
        pval.h2 <- sum(out.mets[2,-1] >= out.mets[2,1], na.rm=TRUE)/(N+1)
        pval.ISA <- sum(out.mets[3,-1] >= out.mets[3,1], na.rm=TRUE)/(N+1)
        pval.mod <- sum(out.mets[4,-1] >= out.mets[4,1], na.rm=TRUE)/(N+1)
        
        
        return(c(nodf=out.mets[1,1], z.nodf=z.nodf, pval.nodf=pval.nodf,
        		h2=out.mets[2,1], z.h2=z.h2, pval.h2=pval.h2,
        		ISA=out.mets[3,1], z.ISA=z.ISA, pval.ISA=pval.ISA,
         		mod=out.mets[4,1], z.mod=z.mod,pval.mod=pval.mod))
      }
    }
  }
  return(rep(NA,12)) 
}
