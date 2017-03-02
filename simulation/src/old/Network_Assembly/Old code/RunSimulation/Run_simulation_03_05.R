dat <-read.csv("C:/Users/Lauren/Desktop/sim_results_02_028.csv")

##  begin simulation

ntree <- 100
ntopo <- 4
nrule <- 5
nmat <- 3

sim.data <- as.data.frame(matrix(0, nrow=1, ncol=11))
names(sim.data) <- names(dat)

for(i in 1:2){
	all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=100, nspecies=30)

	sim.links <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)
	
	sim.links1 <- lapply(sim.links, t)

	sim.links1 <- do.call(rbind, sim.links1)

	sim.links1 <- as.data.frame(sim.links1)

	sim.links1$link.rule <- rownames(sim.links1)

	sim.links1$topo <-rep(c(rep("same.same", nrule*nmat), rep("same.diff", nrule*nmat), rep("diff.diff", nrule*nmat), rep("diff.same", nrule*nmat)), ntree)

	sim.links1$mat <-rep(c("fund", "eco", "samp"), ntree*nrule)

	rownames(sim.links1) <-NULL

	sim.links1$topo <- as.factor(sim.links1$topo)
	sim.links1$link.rule <- as.factor(sim.links1$link.rule)
	sim.links1$mat <- as.factor(sim.links1$mat)
	
	sim.data <- rbind(sim.data, sim.links1)
}
## 


write.table(sim.links, file="sim_results_03_2.csv", sep=",")



