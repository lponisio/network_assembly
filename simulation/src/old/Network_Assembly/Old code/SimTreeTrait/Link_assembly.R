

library(bipartite)
source("/Users/laurenponisio/Documents/R functions/prob_null.R")
source("/Users/laurenponisio/Documents/R functions/break_nest.R")
source("/Users/laurenponisio/Documents/R functions/Mut.adj.R")

link.assembly <- function(p.interacts, model, nsim, fill){
	##p.interacts is a list of probability of interact matrices with dim nplant, npol
	##model can be "abundance", "TM" (traitmatchign), "TM2, "DTM" (descrete traitmatching), "temp" (temporal overlap), "TMTM2" (the two continuous traits), TMTM2TMD (two continuous traits and descrete trait), "ALL" (all the p matrices), "NULL" (forbidden links, temporal overlap and abundance)
 	
		if(model == "abundance"){
			this.interact <- p.interacts$p.abund
		
		} else if (model == "TM"){	
			this.interact <- p.interacts$p.trait
			
		} else if (model == "TM2"){	
			this.interact <- p.interacts$p.trait2
				
		} else if (model == "BAR"){	
			this.interact <- p.interacts$p.trait.bar
			
		} else if (model == "BAR2"){	
			this.interact <- p.interacts$p.trait2.bar
			
		} else if (model == "DTM"){	
			this.interact <- p.interacts$p.discrete.trait
			
		} else if (model == "temp"){	
			this.interact <- p.interacts$p.temp
			
		} else if (model == "TMTM2"){	
			this.interact <- p.interacts$p.trait2 * p.interacts$p.trait
			
		} else if (model == "TMTM2TMD"){	
			this.interact <- p.interacts$p.trait2 * p.interacts$p.trait * p.interacts$p.discrete.trait
		} else if (model == "ALL"){	
			this.interact <- p.interacts$p.trait2 * p.interacts$p.trait * p.interacts$p.discrete.trait * p.interacts$p.abund * p.interacts$p.temp
			
		} else if (model == "NULL"){	
			this.interact <- p.interacts$p.abund * p.interacts$p.temp		
		}				
			
			out.network <- vector(mode="list", length=nsim)
			
			for (i in 1:nsim){
				out.network[[i]] <- matrix(rmultinom(1, fill, prob=as.vector(this.interact)), nrow=nrow(this.interact))
			}
			
	
######calculate nestedness 
		
		nested.network <- break.nest(out.network, method="quasiswap")
		
##calculate degree distribution 

	pol.degree <- vector(mode="list", length=nsim)
	pol.max <- vector(mode="list", length=nsim)
	pol.single <- vector(mode="list", length=nsim)
	
	plant.degree <- vector(mode="list", length=nsim)
	plant.max <- vector(mode="list", length=nsim)
	plant.single <- vector(mode="list", length=nsim)
	
	for (j in 1:nsim){
		pol.degree[[j]] <- apply(out.network[[j]], 2, sum)
		plant.degree[[j]] <- apply(out.network[[j]], 1, sum)
		
		pol.max[[j]] <- max(pol.degree[[j]])
		plant.max[[j]] <- max(plant.degree[[j]])
		
		pol.single[[j]] <- length(pol.degree[[j]][pol.degree[[j]]==1])
		
		plant.single[[j]] <- length(plant.degree[[j]][plant.degree[[j]]==1])
	}
	
	
				
	return(list(networks=out.network, nestedness=nested.network, pol.degree.max=pol.max, plant.degree.max=plant.max, pol.singeltons=pol.single, plant.singeltons=plant.single, pol.degree=pol.degree, plant.degree=plant.degree))
	
}




