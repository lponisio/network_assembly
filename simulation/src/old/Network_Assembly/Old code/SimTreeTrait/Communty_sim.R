library(ape)

library(geiger)

library(distr)


source.comm <- function(ppratio, N ){ ##ppratio is the ratio of plants to pollinators, N is the total number of species
	npol <- ppratio*N
	nplant <- (1-ppratio)*N
 
 ##simulate phylogenies 
	pol.tree <- rcoal(n=npol)
	plant.tree <- rcoal(n=nplant)
 
 ##simulate character evolution
 
 ##discrete
 	###1 is nectar only, 2 is pollen only and 3 is both
	q <- 0.5

	char.pol <- list(matrix(c(-(q+q^2),q,q^2,q,-(q+q^2),q^2,q^2,q^2,-q),3))

	pol.charD <-sim.char(pol.tree, model="discrete", model.matrix=char.pol)
	pol.charD <- pol.charD[,,1] 
	
	char.plant <- list(matrix(c(-q,q,q,-q),2))

	plant.charD <-sim.char(plant.tree, model="discrete", model.matrix=char.plant)
	plant.charD <- plant.charD[,,1] 
	plant.charD <- plant.charD + 1
	
		
	trait.match <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.charD[i]
			for(j in 1:nplant){	
				this.plant <- plant.charD[j]	
				if(this.pol == 3 | this.pol == 2){
					trait.match[j,i] <- 1
				} else if(this.pol == 1 & this.plant==3){
					trait.match[j,i] <- 1
				} else if (this.pol == 1 & this.plant==2){
					trait.match[j,i] <- 0
				}
		}
	}
	
	## define probability matrix
	p.discrete.trait <- trait.match/sum(apply(trait.match,1,sum))
	
##continuous with trait overlap
	##assume proboscus length is log normal, so take exponential 
	
	pol.prob <- exp(rTraitCont(pol.tree, model="BM"))
	plant.corol <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.prob <- cbind((pol.prob - exp(pol.prob)/2), (pol.prob + exp(pol.prob)/2))
	
	range.prob[range.prob < 0] <- 0.01
	
	range.corol <-  cbind((plant.corol - exp(-plant.corol)/2), (plant.corol + exp(-plant.corol)/2))
	
	range.corol[range.corol < 0] <- 0.01
	
	trait.overlap <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.prob[i,]
			for(j in 1:nplant){	
				this.plant <- range.corol[j,]	
				this.overlap <- min(max(this.pol),max(this.plant)) - max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait.overlap[j,i] <- this.overlap
				} else {
					trait.overlap[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait <- trait.overlap/sum(apply(trait.overlap,1,sum))
	
	
	## continuous with trait barrior
	
	trait.bar <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.prob[i]
			for(j in 1:nplant){	
				this.plant <- plant.corol[j]	
				if(this.pol > this.plant){
					trait.bar[j,i] <- 1
				} else {
					trait.bar[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait.bar <- trait.bar/sum(apply(trait.bar,1,sum))
	
	
	## body size
	
	pol.body <- exp(rTraitCont(pol.tree, model="BM"))
	plant.flower <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.body <- cbind((pol.body - exp(-pol.body)/2), (pol.body + exp(-pol.body)/2))
	range.body[range.body < 0] <- 0.01
	
	range.flower <-  cbind((plant.flower - exp(plant.flower)/2), (plant.flower + exp(plant.flower)/2))
	range.flower[range.flower < 0] <- 0.01
	
	body.overlap <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.body[i,]
			for(j in 1:nplant){	
				this.plant <- range.flower[j,]	
				this.overlap <- min(max(this.pol),max(this.plant))- max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					body.overlap[j,i] <- this.overlap
				} else {
					body.overlap[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	p.trait2 <- body.overlap/sum(apply(body.overlap,1,sum))
	
		## continuous with trait barrior
	
	trait2.bar <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.body[i]
			for(j in 1:nplant){	
				this.plant <- plant.flower[j]	
				if(this.pol < this.plant){
					trait2.bar[j,i] <- 1
				} else {
					trait2.bar[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait2.bar <- trait2.bar/sum(apply(trait2.bar,1,sum))
	

########simulate temporal range and define interaction probability 
	
	range.pol <- matrix(runif(2*npol), nrow=npol)
	range.plant <- matrix(runif(2*nplant), nrow=nplant)
	
	temp.overlap <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant[j,]	
				this.overlap <- min(max(this.pol),max(this.plant))- max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					temp.overlap[j,i] <- this.overlap
				} else {
					temp.overlap[j,i] <- 0
				}
		}
	}
	
	
	### interaction probability matrix 
	p.temp <- temp.overlap/sum(apply(temp.overlap,1,sum))
	
##simulate abundances and define interaction probability 

	abund.pol <- ceiling(rlnorm(npol, meanlog=2, sdlog=2))
	abund.plant <- ceiling(rlnorm(nplant, meanlog=2, sdlog=2))
	
	interact.abund <- outer(abund.plant, abund.pol)
	
	## probability of interaction 
	
	p.abund <- interact.abund/sum(apply(interact.abund,1,sum))
	
## equal probability of interating, matrix of 1s
	
	all.1 <- matrix(1, nrow=nplant, ncol=npol)	
	p.equal <- all.1/sum(apply(all.1,1,sum))
	
	return(list(p.abund=p.abund,p.temp=p.temp,p.trait2=p.trait2,p.trait=p.trait,p.discrete.trait=p.discrete.trait, p.trait2.bar=p.trait2.bar, p.trait.bar=p.trait.bar, p.equal=p.equal))
}