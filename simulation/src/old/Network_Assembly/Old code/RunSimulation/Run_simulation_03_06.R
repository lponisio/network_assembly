
##  begin simulation

ntree <- 2
ntopo <- 4
nrule <- 8
nmat <- 5
nloop <- 1
sp <- 30

for(i in 1:nloop){
	all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=100, nspecies=sp)

	sim.links <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)
	
	sim.links1 <- lapply(sim.links, t)

	sim.links1 <- do.call(rbind, sim.links1)

	sim.links1 <- as.data.frame(sim.links1)

	sim.links1$link.rule <- rownames(sim.links1)

	sim.links1$topo <-rep(c(rep("same.same", nrule*nmat), rep("same.diff", nrule*nmat), rep("diff.diff", nrule*nmat), rep("diff.same", nrule*nmat)), ntree)

	sim.links1$mat <-rep(c("fund", "eco", "samp", "samp100", "samp1000"), ntree*nrule)

	rownames(sim.links1) <-NULL

	sim.links1$topo <- as.factor(sim.links1$topo)
	sim.links1$link.rule <- as.factor(sim.links1$link.rule)
	sim.links1$mat <- as.factor(sim.links1$mat)
	
	save(sim.links1, file=sprintf("03_06/%d.RData", i))

}

## 
write.table(sim.data, file="sim_results_02_010.csv", sep=",")