
##  begin simulation

ntree <- 5
ntopo <- 4
nrule <- 5
nmat <- 3

all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=100, nspecies=20)

sim.links4 <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)

## 

sim.links1 <- lapply(sim.links4, t)

sim.links1 <- do.call(rbind, sim.links1)

sim.data <- as.data.frame(sim.links1)

sim.data$link.rule <- rownames(sim.data)

sim.data$topo <-rep(c(rep("same.same", nrule*nmat), rep("same.diff", nrule*nmat), rep("diff.diff", nrule*nmat), rep("diff.same", nrule*nmat)), ntree)

sim.data$mat <-rep(c("fund", "eco", "samp"), ntree*nrule)

rownames(sim.data) <-NULL

sim.data$topo <- as.factor(sim.data$topo)
sim.data$link.rule <- as.factor(sim.data$link.rule)
sim.data$mat <- as.factor(sim.data$mat)

write.table(sim.data, file="sim_results_02_010.csv", sep=",")