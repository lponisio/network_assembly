##  same as `prob.null' but takes `probs' and `fill' as arguments to speed computation
##  and returns quantitative matrix, instead of binary

prob.null2 <- function(probs,fill,nrow,ncol) {
	
	##  draw random interactions
	
	interact <- rmultinom(1,fill,probs)[,1]
	
	return(matrix(interact,nrow=nrow,ncol=ncol))
}

## given inpute matrix (`data') returns desired summary statistics

calc.metric <- function(data) {
	
	data <- as.matrix(empty(data))
	
	mets <- networklevel(data, index=c("weighted NODF", "ISA", "H2"), H2_integer=TRUE)
	print("calculated mets")

	return(mets)
}


## function to simulate 1 null, and calculate statistics on it
null.stat <- function(prob,fill,nr,nc) {
	print("null_sim")
	
	## simulate null web, previous versions used r2dtable or vaznull.fast
		sim.web <- prob.null2(prob,fill,nr,nc)
		
	####following line is only necessary when using a null model that does not constraint marginals, i.e. prob.null2
		while(any(dim(sim.web) < 5) ) sim.web <- prob.null2(prob,fill,nr,nc) 
	
	return(calc.metric(sim.web))
}


##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)

network.metrics <- function (data, N) {
    
	print(data)
	 ## check that matrix is proper format (no empty row/col and no NAs)
    if(sum(data, na.rm=TRUE) != 0 & !any(is.na(data))) {
		print("1")
		
	##round decimals to integers for calculation of statistics
		if(any(data[data != 0] < 1)) {
    		data <- round(data*10) #(1/min(data[data != 0]))
    	}
    	print("2")
    ## drop empty rows and columns
    	data <- as.matrix(empty(data))
		print("3")
		
	## check to make sure the rounded,emptied matrix is large enough to calculate statistics on
		if(dim(data)[1] > 5 & dim(data)[2] > 5) {
		print("4")
    	
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
			#pval.mod <- sum(out.mets[5,-1] >= out.mets[5,1], na.rm=TRUE)/(N+1)
		
			return(c(nofd = out.mets[1,1], pval.nodf=pval.nodf,
					ISA=out.mets[2,1], pval.ISA=pval.ISA, 
					h2=out.mets[3,1], pval.h2= pval.h2 )) 
					#mod = out.mets[5,1], pval.mod=pval.mod,
					#num.mod = out.mets[4,1]))
		} else {
	
			return(rep(NA,6)) 
		}

	} else {
	
		return(rep(NA,6)) ##change back to 9 if considering modularity and number of modules
	
	}
}
