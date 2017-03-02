library(ape)

library(geiger)

## tree and trait can be "same" or "diff"  
## if tree == "same" tree of plants and animals is created separatly
## if tree == "diff" tree of plants and animals is identical 
## if trait == "same" , trait.type must be "cont" 
## the character evolution is simulated for plants and then animal traits are drawn form the plant traits
## if trait =="diff", trait.type can be "cont" or "discrete"
## trait evolution is simulated separtly for animals and plants

phlyo.comm <- function(tree, traits, trait.type, N){
	
	
	if(tree=="diff"){
		
		plant.tree <- rcoal(n=N)
		animal.tree <- rcoal(n=N)
		
	} else if (tree == "same"){
		
	plant.tree <- animal.tree  <- rcoal(n=N)
	
	}
	
	if (traits == "diff"){
		
		if(trait.type="cont"){
	
			animal.trait <- exp(rTraitCont(animal.tree, model="BM"))
			plant.trait <- exp(rTraitCont(plant.tree, model="BM"))
	
		} else if (trait.type="discrete"){
		
			###1 is nectar only, 2 is pollen only and 3 is both
 			## p(transitioning between states) is equal any way
 			##plants can take on 1,3 pollinators 1,2,3
	
			q <- 0.5
	
			## pollinator character evolution 
	
			animal.state <- list(matrix(c(-(q+q^2),q,q^2,q,-(q+q^2),q^2,q^2,q^2,-q),3))

			animal.trait <-sim.char(animal.tree, model="discrete", model.matrix=animal.state)
			animal.trait <- animal.trait[,,1] 
	
		## plant character evolution 
	
			plant.state <- list(matrix(c(-q,q,q,-q),2))

			plant.trait <-sim.char(plant.tree, model="discrete", model.matrix=char.plant)
			plant.trait <- plant.charD[,,1] + 1

		
		}
		
	
	} else if (traits == "same"){

		plant.trait <- exp(rTraitCont(plant.tree, model="BM"))
		animal.trait <- sample(plant.trait, N)
		
	}
	
	return(list(plant.trait=plant.trait, animal.trait=animal.trait))
}

