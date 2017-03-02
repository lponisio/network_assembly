library(ape)
library(geiger)
library(distr)

source("/Users/laurenponisio/Documents/Network Assembly /metrics_networks_sim.R")

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

## ranges are equally limited (wide)

	
	trait.wide <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol.wide[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant.wide[j,]	
				this.overlap <- min(max(this.pol),max(this.plant)) - max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait.wide[j,i] <- this.overlap
				} else {
					trait.wide[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait.wide <- trait.wide/sum(apply(trait.wide,1,sum))
	
	
## ranges are equally limited (narrow)

		
	trait.narrow <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol.narrow[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant.narrow[j,]	
				this.overlap <- min(max(this.pol),max(this.plant)) - max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait.narrow[j,i] <- this.overlap
				} else {
					trait.narrow[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait.narrow <- trait.narrow/sum(apply(trait.narrow,1,sum))
	

##larger trait size limits range for plant, increases from pol 
	
	trait.comp <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol.wide[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant.narrow[j,]	
				this.overlap <- min(max(this.pol),max(this.plant)) - max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait.comp[j,i] <- this.overlap
				} else {
					trait.comp[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait.comp <- trait.comp/sum(apply(trait.comp,1,sum))
	
	
	## continuous with trait barrior
	
	trait.bar <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.trait[i]
			for(j in 1:nplant){	
				this.plant <- plant.trait[j]	
				if(this.pol > this.plant){
					trait.bar[j,i] <- 1
				} else {
					trait.bar[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait.bar <- trait.bar/sum(apply(trait.bar,1,sum))
	
	## equal probability of interating, matrix of 1s
	
	neutral <- matrix(1, nrow=nplant, ncol=npol)	
	p.neutral <- neutral/sum(apply(neutral,1,sum))
	
##simulate abundances and define interaction probability 

	abund.pol <- ceiling(rlnorm(npol, meanlog=2, sdlog=2))
	abund.plant <- ceiling(rlnorm(nplant, meanlog=2, sdlog=2))
	
	interact.abund <- outer(abund.plant, abund.pol)
	
	## probability of interaction 
	
	p.abund <- interact.abund/sum(apply(interact.abund,1,sum))
	

### group together funmamental p matrices 

	 fundamental <- list(
	 	trait.comp=trait.comp, 
	 	trait.bar=trait.bar,
	 	neutral=neutral, 
	 	trait.narrow= trait.narrow,
	 	trait.wide=trait.wide)
	 	 
	 p.fundamental <- list(
	 	p.trait.comp = p.trait.comp, 
	 	p.trait.bar = p.trait.bar,
	 	p.neutral = p.neutral, 
	 	p.trait.narrow = p.trait.narrow,
	 	p.trait.wide = p.trait.wide) 	 
	
	p.realized <- lapply(p.fundamental, FUN = function(X){X * p.abund})
	
	realized <- vector(mode="list", length=length(p.realized))
			
	#for(j in 1:length(p.realized)){
		#for (i in 1:nsim){
				#realized[[j]][[i]] <- matrix(rmultinom(1, nint, prob=as.vector(p.realized[[j]])), nrow=nrow(p.realized[[j]]))
		#}
	#}

	nested.fundamental <- metric.net(fundamental, null="quasiswap", func=nestednodf, nsim=100) 	
	#h2.fundamental <- metric.net(fundamental, null="quasiswap", func=h2, nsim=100) 	
	
	return(nested.fundamental)
	
} ##close function