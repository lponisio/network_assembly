
library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)
library(sna)

##creat graph from matrix

mut.adj <- function(x) { ##matrix x
		nr <- dim(x)[1]
		nc <- dim(x)[2]
		
		to.fill <- matrix(0,ncol = nc + nr, nrow = nc + nr)
		
		to.fill[1:nr,(nr+1):(nc+nr)] <- x
		
		adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
		
		return(adj.mat)
}

###simulate trees

sim.phylo <- function(nsim, mu, lambda, age, nspecies){
	
	tree.1 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	tree.2 <- sim.bd.taxa.age(n = nspecies, numbsim = nsim, lambda = lambda, mu = mu, frac = 1, age=age, mrca = TRUE)
	
	return(list(tree.1, tree.2))
}


##simulate character evolution
character.evo <- function(tree1, tree2){
	
	tree1.trait1 <- exp(rTraitCont(tree1, model="BM"))
	tree1.trait2 <- exp(rTraitCont(tree1, model="BM"))
	
	tree2.trait1 <- exp(rTraitCont(tree2, model="BM"))
	
	sample.trait <- sample(log(tree1.trait1), size = Ntip(tree2), replace=TRUE)
	
	tree2.trait2 <- as.vector(exp(rmvnorm(1,sample.trait,sigma=0.05*vcv(tree2),method="chol")))
	
	same.same <- list(tree1.trait1, tree1.trait1)
	same.diff <- list(tree1.trait1, tree1.trait2)		
	diff.diff <- list(tree1.trait1, tree2.trait1)
	diff.same <- list(tree1.trait1, tree2.trait2)	
	
	
	return(list(same.same = same.same, same.diff = same.diff, diff.diff = diff.diff, diff.same = diff.same))
}


##bring everything together

master.fun <- function(tree1,tree2) {
	
	these.trait <- character.evo(tree1, tree2)
	
	same.same.stat <- link.rule(these.trait$same.same, nsim=100)
	same.diff.stat <- link.rule(these.trait$same.diff, nsim=100)
	diff.same.stat <- link.rule(these.trait$diff.same, nsim=100)
	diff.diff.stat <- link.rule(these.trait$diff.diff, nsim=100)

	return(list(same.same.stat,same.diff.stat,diff.same.stat,diff.diff.stat))
}


###simulate linkages



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

####calculate statistics


network.metrics <- function(data){
		
	if(sum(data != 0)){
		
		if(any(data < 1)){
		
			data <- round(data*(100/min(data[data !=0])))
		}
		
	
	
		data <- empty(data, count=FALSE)
		
		nulls <- vaznull(N=1, web=data)
				
		
		nets <- c(list(data=data), nulls=nulls)
	
		out.mets <- sapply(nets, function(X){

			mets <- networklevel(X, index=c("weighted NODF", "ISA", "H2"), H2_integer=FALSE)
	
			graph <- mut.adj(X)
			#plot(graph)
			weights <- as.vector(X)
			weights <- weights[weights != 0]
	
			wtc <- walktrap.community(graph, weights=weights, steps=1000)
			memb <- community.to.membership(graph, wtc$merges, steps= (nrow(wtc$merges) -1) )
			mod.met <- modularity(graph, memb$membership, weights=weights)
			num.mods <- length(unique(memb$membership))
				
			return(c(mets, num.mods=num.mods, mod.met=mod.met))

		})## end out.mets apply function
	
		pval.nodf <- sum(out.mets[1,-1] >= out.mets[1,"data"])/length(out.mets)
		pval.ISA <- sum(out.mets[2,-1] >= out.mets[2,"data"])/length(out.mets)
		pval.h2 <- sum(out.mets[3,-1] >= out.mets[3,"data"])/length(out.mets)
		pval.mod <- sum(out.mets[5,-1] >= out.mets[5,"data"])/length(out.mets)
	
		return(c(nofd = out.mets[1,"data"], pval.nodf=pval.nodf,
				ISA=out.mets[2,"data"], pval.ISA=pval.ISA, 
				h2=out.mets[3,"data"], pval.h2= pval.h2, 
				mod= out.mets[5,"data"],  pval.mod=pval.mod, num.mod= out.mets[4,"data"]))
			 
	} else{ 
		
	return(c(rep(NA, 9)))
	}
}



####

ntree <- 1
ntopo <- 4
nrule <- 5

all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=100, nspecies=10)

sim.links <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)

sim.stats <-rapply(sim.links, f=network.metrics, how="replace")

sim.data <- t(matrix(rapply(sim.stats, c),nrow=9))

sim.data<- cbind(rep(c("trait.comp", "trait.bar", "neutral", "trait.narrow", "trait.wide"), ntopo*ntree), sim.data)
	
	
sim.data <- cbind(rep(c(rep("same.same", nrule), rep("same.diff", nrule), rep("diff.diff", nrule), rep("diff.same", nrule)), ntree), sim.data)
	
	
colnames(sim.data) <- c( "topo", "linkr", "nofd", "pval.nodf", "ISA", "pval.ISA", "h2", "pval.h2", "mod", "pval.mod", "num.mod")

save.image("./sim_out_test.Rdata")


