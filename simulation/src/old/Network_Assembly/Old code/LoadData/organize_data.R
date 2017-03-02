

organize.data <- function(data, ntree, ntopo, nrule){
	
	out.data <- matrix(NA, nrow=ntree*ntopo*nrule, ncol=9)

	
	for(i in 1:ntree){
		for(j in 1:ntopo){
			for(l in 1:nrule){
				for(z in 1: (ntree*ntopo*nrule)){
					out.data[z] <- data[[i]][[j]][[l]]
					
				}
			}
		}
	}
	
		
	out.data <- cbind(rep(c("trait.comp", "trait.bar", "neutral", "trait.narrow", "trait.wide"), ntopo*ntree), out.data)
	
	
	out.data <- cbind(rep(c(rep("same.same", nrule), rep("same.diff", nrule), rep("diff.diff", nrule), rep("diff.same", nrule)), ntree), out.data)
	
	
	colnames(out.data) <- c( "topo", "linkr", "nofd", "pval.nodf", "ISA", "pval.ISA", "h2", "pval.h2", "mod", "pval.mod", "num.mod")

	
	return(out.data)
}