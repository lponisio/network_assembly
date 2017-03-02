library(ape)

library(geiger)

library(distr)


source.comm <- function(ppratio, N, nsim, nint ){ ##ppratio is the ratio of plants to pollinators, N is the total number of species
	npol <- ppratio*N
	nplant <- (1-ppratio)*N
 
 ##simulate phylogenies 
	pol.tree <- rcoal(n=npol)
	plant.tree <- rcoal(n=nplant)
 
 
 ###define interaction matrices based on traits
 
 ## simulate discrete character evolution
 	###1 is nectar only, 2 is pollen only and 3 is both
 	## p(transitioning between states) is equal any way
 	##plants can take on 1,3 pollinators 1,2,3
	
	q <- 0.5
	
	## pollinator character evolution 
	
	char.pol <- list(matrix(c(-(q+q^2),q,q^2,q,-(q+q^2),q^2,q^2,q^2,-q),3))

	pol.charD <-sim.char(pol.tree, model="discrete", model.matrix=char.pol)
	pol.charD <- pol.charD[,,1] 
	
	## plant character evolution 
	
	char.plant <- list(matrix(c(-q,q,q,-q),2))

	plant.charD <-sim.char(plant.tree, model="discrete", model.matrix=char.plant)
	plant.charD <- plant.charD[,,1] 
	plant.charD <- plant.charD + 1
	
		
	discrete.trait <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.charD[i]
			for(j in 1:nplant){	
				this.plant <- plant.charD[j]	
				if(this.pol == 3 | this.pol == 2){
					discrete.trait[j,i] <- 1
				} else if(this.pol == 1 & this.plant==3){
					discrete.trait[j,i] <- 1
				} else if (this.pol == 1 & this.plant==2){
					discrete.trait[j,i] <- 0
				}
		}
	}
	
	## define probability matrix
	
	p.discrete.trait <- discrete.trait/sum(apply(discrete.trait,1,sum))
	
##continuous with trait overlap

## ranges are equally limited (wide)

	pol.trait <- exp(rTraitCont(pol.tree, model="BM"))
	plant.trait <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.pol <- cbind((pol.trait - exp(pol.trait)/2), (pol.trait + exp(pol.trait)/2))
	
	range.pol[range.pol < 0] <- 0.01
	
	range.plant <-  cbind((plant.trait - exp(plant.trait)/2), (plant.trait + exp(plant.trait)/2))
	
	range.plant[range.plant < 0] <- 0.01
	
	trait.wide <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant[j,]	
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

	pol.trait <- exp(rTraitCont(pol.tree, model="BM"))
	plant.trait <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.pol <- cbind((pol.trait - exp(-pol.trait)/2), (pol.trait + exp(-pol.trait)/2))
	
	range.pol[range.pol < 0] <- 0.01
	
	range.plant <-  cbind((plant.trait - exp(-plant.trait)/2), (plant.trait + exp(-plant.trait)/2))
	
	range.plant[range.plant < 0] <- 0.01
	
	trait.narrow <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant[j,]	
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
	

##larger tait size limits range for plant, increses from pol 

	##assume proboscus length is log normal, so take exponential 
	
	pol.prob <- exp(rTraitCont(pol.tree, model="BM"))
	plant.corol <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.prob <- cbind((pol.prob - exp(pol.prob)/2), (pol.prob + exp(pol.prob)/2))
	
	range.prob[range.prob < 0] <- 0.01
	
	range.corol <-  cbind((plant.corol - exp(-plant.corol)/2), (plant.corol + exp(-plant.corol)/2))
	
	range.corol[range.corol < 0] <- 0.01
	
	trait1 <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.prob[i,]
			for(j in 1:nplant){	
				this.plant <- range.corol[j,]	
				this.overlap <- min(max(this.pol),max(this.plant)) - max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait1[j,i] <- this.overlap
				} else {
					trait1[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait1 <- trait1/sum(apply(trait1,1,sum))
	
	
	## continuous with trait barrior
	
	pol.prob <- exp(rTraitCont(pol.tree, model="BM"))
	plant.corol <- exp(rTraitCont(plant.tree, model="BM"))
	
	trait1.bar <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- pol.prob[i]
			for(j in 1:nplant){	
				this.plant <- plant.corol[j]	
				if(this.pol > this.plant){
					trait1.bar[j,i] <- 1
				} else {
					trait1.bar[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	
	p.trait1.bar <- trait1.bar/sum(apply(trait1.bar,1,sum))
	
	
	##larger trait size limits range for pol, increased for plant 	
	pol.body <- exp(rTraitCont(pol.tree, model="BM"))
	plant.flower <- exp(rTraitCont(plant.tree, model="BM"))
	
	range.body <- cbind((pol.body - exp(-pol.body)/2), (pol.body + exp(-pol.body)/2))
	range.body[range.body < 0] <- 0.01
	
	range.flower <-  cbind((plant.flower - exp(plant.flower)/2), (plant.flower + exp(plant.flower)/2))
	range.flower[range.flower < 0] <- 0.01
	
	trait2 <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.body[i,]
			for(j in 1:nplant){	
				this.plant <- range.flower[j,]	
				this.overlap <- min(max(this.pol),max(this.plant))- max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					trait2[j,i] <- this.overlap
				} else {
					trait2[j,i] <- 0
				}
		}
	}
	
	### interaction probability matrix 
	p.trait2 <- trait2/sum(apply(trait2,1,sum))
	
		## continuous with trait barrior
	
	pol.body <- exp(rTraitCont(pol.tree, model="BM"))
	plant.flower <- exp(rTraitCont(plant.tree, model="BM"))
	
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
	
	temp <- matrix(NA, nrow=nplant, ncol=npol)
	
	for(i in 1:npol){
			this.pol <- range.pol[i,]
			for(j in 1:nplant){	
				this.plant <- range.plant[j,]	
				this.overlap <- min(max(this.pol),max(this.plant))- max(min(this.pol),min(this.plant)) 
				if(this.overlap > 0){
					temp[j,i] <- this.overlap
				} else {
					temp[j,i] <- 0
				}
		}
	}
	
	
	### interaction probability matrix 
	p.temp <- temp/sum(apply(temp,1,sum))
	
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

	 fundamental <- list(temp=temp,
	 	trait2=trait2,
	 	trait1=trait1,
	 	discrete.trait=discrete.trait,
	 	trait2.bar=trait2.bar, 
	 	trait1.bar=trait1.bar,
	 	neutral=neutral, 
	 	trait.narrow= trait.narrow,
	 	trait.wide=trait.wide)
	 	 
	 p.fundamental <- list(p.temp=p.temp,
	 	p.trait2=p.trait2,
	 	p.trait1=p.trait1,
	 	p.discrete.trait=p.discrete.trait,
	 	p.trait2.bar=p.trait2.bar, 
	 	p.trait1.bar=p.trait1.bar,
	 	p.neutral=p.neutral, 
	 	p.trait.narrow= p.trait.narrow,
	 	p.trait.wide=p.trait.wide) 	 
	
	p.realized <- lapply(p.fundamental, FUN = function(X){X * p.abund})
	
	realized <- vector(mode="list", length=length(p.realized))
			
	for(j in 1:length(p.realized)){
		for (i in 1:nsim){
				realized[[j]][[i]] <- matrix(rmultinom(1, nint, prob=as.vector(p.realized[[j]])), nrow=nrow(p.realized[[j]]))
		}
	}

	#nested.realized <- metric.net(realized, null="quasiswap", func=nestednodf) 	
	motifs.realized <- rapply(realized, f=bmotifs, how="replace")
	
	return(motifs.realized)
	
} ##close function