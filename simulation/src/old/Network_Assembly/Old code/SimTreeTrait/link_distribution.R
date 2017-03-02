##simulate linkage distribution

	##truncated power law 
	#t.pl <- function(k,y,kc){
		#p <- (k^-y)*exp(-k/kc)
		#return(p)
	#}

	#pol.y <- plant.y <- 1
	#pol.kc <- nplant*0.8
	#plant.kc <- npol*0.8

	#degree.dist.pol <- DiscreteDistribution(1:nplant, prob=(1/sum(t.pl(y=pol.y, k=1:nplant,kc=pol.kc)))*t.pl(y=pol.y, k=1:nplant,kc=pol.kc))
	
	
	#degree.dist.plant <- DiscreteDistribution(1:nplant, prob=(1/sum(t.pl(y=plant.y, k=1:npol,kc=plant.kc)))*t.pl(y=plant.y, k=1:npol,kc=plant.kc))


	#pol.link <- degree.dist.pol@r(npol)
	#plant.link <- degree.dist.plant@r(nplant)