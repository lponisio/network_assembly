
library(igraph)

library(bipartite)



source("/Users/laurenponisio/Documents/Hedgerow Networks/Code/Functions/Mut.adj.R")

network.metrics <- function(data){
		
	if(sum(data != 0)){
		
		if(any(data < 1)){
		
			data <- round(data*(100/min(data[data !=0])))
		}
	
		data <- empty(data, count=FALSE)
		
		nulls <- vaznull(N=100, web=data)
			
		nets <- c(list(data=data), nulls=nulls)
	
		out.mets <- sapply(nets, function(X){
			
			mets <- networklevel(X, index=c("weighted NODF", "ISA", "H2"), H2_integer=FALSE)
	
			graph <- mut.adj(X)
			#plot(graph)
			weights <- as.vector(X)
			weights <- weights[weights != 0]
	
			wtc <- walktrap.community(graph, weights=weights, steps=1000)
			memb <- community.to.membership(graph, wtc$merges, steps= (nrow(wtc$merges) -1) )
			mod.met <- modularity(graph, memb$membership, weights=weights)
			num.mods <- length(unique(memb$membership))
				
			return(c(mets, num.mods=num.mods, mod.met=mod.met))

		})## end out.mets apply function
	
		pval.nodf <- sum(out.mets[1,-1] >= out.mets[1,"data"])/ncol(out.mets)
		pval.ISA <- sum(out.mets[2,-1] >= out.mets[2,"data"])/ncol(out.mets)
		pval.h2 <- sum(out.mets[3,-1] >= out.mets[3,"data"])/ncol(out.mets)
		pval.mod <- sum(out.mets[5,-1] >= out.mets[5,"data"])/ncol(out.mets)
	
		return(c(nofd = out.mets[1,"data"], pval.nodf=pval.nodf,
				ISA=out.mets[2,"data"], pval.ISA=pval.ISA, 
				h2=out.mets[3,"data"], pval.h2= pval.h2, 
				mod= out.mets[5,"data"],  pval.mod=pval.mod, num.mod= out.mets[4,"data"]))
			 
	} else{ 
		
		return(rep(NA, 9))
	}
}

