
interact.sample <- function(abundance, trait.overlap, npol, nplant){
	
	trait <- matrix(trait.overlap, nrow=npol)
	
	if(sum(trait) != 0){
	
		trait.abund <- trait*abundance
	
		fill <- round(sum(trait.abund))
	
		p.trait.abund <- trait.abund/sum(apply(trait.abund,1,sum))
	
		trait.samp <- matrix(rmultinom(1, fill, prob=as.vector(p.trait.abund)), nrow=nrow(p.trait.abund))
	
		return(list(trait, trait.abund, trait.samp))
	
	} else{
		
		x <- matrix(0, nrow=npol, ncol=nplant)
		
		return(list(x,x,x))
	}	
}