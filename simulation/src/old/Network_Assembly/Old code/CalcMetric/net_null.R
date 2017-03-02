## given inpute matrix (`data') returns desired summary statistics

calc.metric <- function(data) {
	
	data <- as.matrix(empty(data))
	
	mets <- networklevel(data, index=c("weighted NODF", "ISA", "H2"), H2_integer=TRUE)
	print("calculated mets")
	graph <- mut.adj(data)
	print("converted to graph")
#	plot(graph)
	weights <- as.vector(data)
	weights <- weights[weights != 0]
	
	wtc <- walktrap.community(graph, weights=weights, steps=100)
	memb <- community.to.membership(graph, wtc$merges, steps= (nrow(wtc$merges) -1))
	mod.met <- modularity(graph, memb$membership, weights=weights)
	num.mods <- length(unique(memb$membership))
	
	print("calculated modularity")
	return(c(mets, num.mods=num.mods, mod.met=mod.met))
}


## function to simulate 1 null, and calculate statistics on it
null.stat <- function(prob,fill,nr,nc) {
	print("null_sim")
	
	## simulate null web, previous versions used r2dtable or vaznull.fast
	sim.web <- prob.null2(prob,fill,nr,nc)
	
	return(calc.metric(sim.web))
}


##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)

network.metrics <- function (data, N) {
    ## check that matrix is proper format (no empty row/col)
   
    if(sum(data) != 0) {
    	
    	if(any(data[data != 0] < 1)) {
    		data <- round(data*10) #(1/min(data[data != 0]))
    	}
    	data <- as.matrix(empty(data))
    	
    	print(data)
    	
    	##  produce N sets of statistics form N different nulls
    	##  should return matrix that is 5 row by N col
		
		row.prob <- rowSums(data)
		col.prob <- colSums(data)
		data.prob <- expand.grid(row.prob,col.prob)
		data.prob <- data.prob[,1]*data.prob[,2]
		
		data.fill <- sum(data)
		
		null.stat <- replicate(N, null.stat(data.prob,data.fill,nrow(data),ncol(data)), simplify = TRUE)
		print("true_met")
		true.stat <- calc.metric(data)
		
		##  combind into one matrix from which p-values will be calculated
		out.mets <- cbind(true.stat,null.stat)
		
		##  compute p-vals
		pval.nodf <- sum(out.mets[1,-1] >= out.mets[1,1], na.rm=TRUE)/(N+1)
		pval.ISA <- sum(out.mets[2,-1] >= out.mets[2,1], na.rm=TRUE)/(N+1)
		pval.h2 <- sum(out.mets[3,-1] >= out.mets[3,1], na.rm=TRUE)/(N+1)
		pval.mod <- sum(out.mets[5,-1] >= out.mets[5,1], na.rm=TRUE)/(N+1)
		
		return(c(nofd = out.mets[1,1], pval.nodf=pval.nodf,
				ISA=out.mets[2,1], pval.ISA=pval.ISA, 
				h2=out.mets[3,1], pval.h2= pval.h2, 
				mod = out.mets[5,1], pval.mod=pval.mod,
				num.mod = out.mets[4,1]))
	
	} else {
	
		return(rep(NA,9))
	
	}
}



##creat graph from matrix, needed to use modularity optimization function

mut.adj <- function(x) { ##matrix x
		nr <- dim(x)[1]
		nc <- dim(x)[2]
		
		to.fill <- matrix(0,ncol = nc + nr, nrow = nc + nr)
		
		to.fill[1:nr,(nr+1):(nc+nr)] <- x
		
		adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
		
		return(adj.mat)
}


##  same as `prob.null' but takes `probs' and `fill' as arguments
##  and returns quantitative matrix, instead of binary
prob.null2 <- function(probs,fill,nrow,ncol) {
	
	##  draw random interactions
	interact <- rmultinom(1,fill,probs)[,1]
	
	return(matrix(interact,nrow=nrow,ncol=ncol))
}


#####################  old null simulator  #####################


## fast function to simulate 1 null matrix
vaznull.fast <- function(web) {
	rs.p <- rowSums(web)/sum(web)
	cs.p <- colSums(web)/sum(web)
	P <- P1 <- tcrossprod(rs.p, cs.p)
	finalmat <- matrix(0, nrow(web), ncol(web))
	
	n.int.finalmat <- 0
	i <- 1
	while (n.int.finalmat < sum(dim(web))) {
		if(i == 10000) print("watch out!!!!")
		sel <- sample(1:length(web), 1, prob = P, replace = TRUE)
		selc <- floor((sel - 1)/(dim(web)[1])) + 1
		selr <- ((sel - 1)%%dim(web)[1]) + 1
		if (sum(finalmat[, selc]) == 0 | sum(finalmat[selr,]) == 0) {
			finalmat[sel] <- 1
			P[sel] <- 0
		}
		
		n.int.finalmat <- sum(rowSums(finalmat) > 0) + sum(colSums(finalmat) > 0)
        i <- 1+i
	}
	
	conn.remain <- sum(web > 0) - sum(finalmat > 0)
	if (conn.remain > 0) {
		add <- sample(which(finalmat == 0), conn.remain, prob = P1[finalmat == 0])
		finalmat[add] <- 1
	}
	
	int.remain <- sum(web) - sum(finalmat)
	if (int.remain > 0) {
		add <- sample(which(finalmat > 0), int.remain, prob = P1[finalmat > 0], replace = TRUE)
		finalmat[as.numeric(names(table(add)))] <- finalmat[as.numeric(names(table(add)))] + table(add)
	}
	
	return(finalmat)
}