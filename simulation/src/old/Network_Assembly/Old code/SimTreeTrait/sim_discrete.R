discrete.evo <- function(tree1, tree2){
	###1 is nectar only, 2 is pollen only and 3 is both
 	## p(transitioning between states) is equal any way
 	##plants can take on 1,3 pollinators 1,2,3
	
	q <- 0.5
	
	## pollinator character evolution 
	
	animal.state <- list(matrix(c(-(q+q^2),q,q^2,q,-(q+q^2),q^2,q^2,q^2,-q),3))

	tree1.d <-sim.char(tree1, model="discrete", model.matrix=animal.state)[,,1] 
	
	## plant character evolution 
	
	plant.state <- list(matrix(c(-q,q,q,-q),2))

	tree2.d <- sim.char(tree2, model="discrete", model.matrix=char.plant)[,,1] + 1
	
	discrete.trait <- matrix(NA, nrow=Ntip(tree1), ncol=Ntip(tree2))
	
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
	
	
}	

