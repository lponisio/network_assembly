library(ape)
library(geiger)
#library(distr)



link.rule <- function(traits, nsim){ ##ppratio is the ratio of plants to pollinators, N is the total number of species

	npol <- nplant <- length(traits[[1]])
	
	pol.trait <- traits[[1]] 
	plant.trait <- traits[[2]]
	
##continuous with trait overlap

	##wide 
	range.pol.wide <- cbind((pol.trait - exp(pol.trait)/2), (pol.trait + exp(pol.trait)/2))
	range.pol.wide[range.pol.wide < 0] <- 0.01
	
	range.plant.wide <-  cbind((plant.trait - exp(plant.trait)/2), (plant.trait + exp(plant.trait)/2))
	range.plant.wide[range.plant.wide < 0] <- 0.01

	#narrow
	range.pol.narrow <- cbind((pol.trait - exp(-pol.trait)/2), (pol.trait + exp(-pol.trait)/2))
	range.pol.narrow[range.pol.narrow < 0] <- 0.01
	
	range.plant.narrow <-  cbind((plant.trait - exp(-plant.trait)/2), (plant.trait + exp(-plant.trait)/2))
	range.plant.narrow[range.plant.narrow < 0] <- 0.01

	pa.comb <- expand.grid(1:npol , 1:nplant)
	
	## narrow traits 
	pol.narrow.min <- range.pol.narrow[pa.comb[,1],1]
	pol.narrow.max <- range.pol.narrow[pa.comb[,1],2]
	
	plant.narrow.min <- range.plant.narrow[pa.comb[,2],1]
	plant.narrow.max <- range.plant.narrow[pa.comb[,2],2]

	## wide traits
	
	pol.wide.min <- range.pol.wide[pa.comb[,1],1]
	pol.wide.max <- range.pol.wide[pa.comb[,1],2]
	
	plant.wide.min <- range.plant.wide[pa.comb[,2],1]
	plant.wide.max <- range.plant.wide[pa.comb[,2],2]

	
	##narrow overlap
	
	lower.lim.narrow <- ifelse(pol.narrow.min > plant.narrow.min, pol.narrow.min, plant.narrow.min)
	upper.lim.narrow <- ifelse(pol.narrow.max < plant.narrow.max, pol.narrow.max, plant.narrow.max)
	
	trait.overlap.narrow <- upper.lim.narrow - lower.lim.narrow 
	trait.overlap.narrow[trait.overlap.narrow < 0] <- 0
	
	trait.narrow <- matrix(trait.overlap.narrow, nrow=npol)
	
	p.trait.narrow <- trait.narrow/sum(apply(trait.narrow,1,sum))
	
	#wide overlap
	
	lower.lim.wide <- ifelse(pol.wide.min > plant.wide.min, pol.wide.min, plant.wide.min)
	upper.lim.wide <- ifelse(pol.wide.max < plant.wide.max, pol.wide.max, plant.wide.max)
	
	trait.overlap.wide <- upper.lim.wide - lower.lim.wide 
	trait.overlap.wide[trait.overlap.wide < 0] <- 0
	
	trait.wide <- matrix(trait.overlap.wide, nrow=npol)
	
	p.trait.wide <- trait.wide/sum(apply(trait.wide,1,sum))
	
	##one wide, one narrow
	
	lower.lim.comp <- ifelse(pol.wide.min > plant.narrow.min, pol.wide.min, plant.narrow.min)
	upper.lim.comp <- ifelse(pol.wide.max < plant.narrow.max, pol.wide.max, plant.narrow.max)
	
	trait.overlap.comp <- upper.lim.comp - lower.lim.comp 
	trait.overlap.comp[trait.overlap.comp < 0] <- 0
	
	trait.comp <- matrix(trait.overlap.comp, nrow=npol)

	p.trait.comp <- trait.comp/sum(apply(trait.comp,1,sum))
	
	
	## continuous with trait barrior
	
	plant.trait.combin <- plant.trait[pa.comb[,2]]
	pol.trait.combin <- plant.trait[pa.comb[,1]]
	
	barrior <- ifelse(pol.trait.combin > plant.trait.combin, 1, 0)
	
	trait.bar <- matrix(barrior, nrow=npol)
		
	p.trait.bar <- trait.bar/sum(apply(trait.bar,1,sum))
	
	## equal probability of interating, matrix of 1s
	
	neutral <- matrix(1, nrow=nplant, ncol=npol)	
	p.neutral <- neutral/sum(apply(neutral,1,sum))
	
##simulate abundances and define interaction probability 

	#abund.pol <- ceiling(rlnorm(npol, meanlog=2, sdlog=2))
	#abund.plant <- ceiling(rlnorm(nplant, meanlog=2, sdlog=2))
	
	#interact.abund <- outer(abund.plant, abund.pol)
	
	## probability of interaction 
	
	#p.abund <- interact.abund/sum(apply(interact.abund,1,sum))
	

### group together funmamental p matrices 

	 fundamental <- list(
	 	trait.comp=trait.comp, 
	 	trait.bar=trait.bar,
	 	neutral=neutral, 
	 	trait.narrow= trait.narrow,
	 	trait.wide=trait.wide)
	 	 
	 #p.fundamental <- list(
	 	#p.trait.comp = p.trait.comp, 
	 	#p.trait.bar = p.trait.bar,
	 	#p.neutral = p.neutral, 
	 	#p.trait.narrow = p.trait.narrow,
	 	#p.trait.wide = p.trait.wide) 	 
	
	#p.realized <- lapply(p.fundamental, FUN = function(X){X * p.abund})
	
	#realized <- vector(mode="list", length=length(p.realized))
			
	#for(j in 1:length(p.realized)){
		#for (i in 1:nsim){
				#realized[[j]][[i]] <- matrix(rmultinom(1, nint, prob=as.vector(p.realized[[j]])), nrow=nrow(p.realized[[j]]))
		#}
	#}


	
	return(fundamental)
	
} ##close function



