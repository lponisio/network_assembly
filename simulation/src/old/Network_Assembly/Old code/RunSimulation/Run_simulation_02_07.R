
##  begin simulation

ntree <- 100
ntopo <- 4
nrule <- 4

all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=100, nspecies=30)

sim.links2 <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)

## 

sim.links1 <- lapply(sim.links2, t)

sim.links1 <- do.call(rbind, sim.links1)

sim.data <- as.data.frame(sim.links1)

sim.data$link.rule <- rownames(sim.data)

sim.data$topo <-rep(c(rep("same.same", nrule), rep("same.diff", nrule), rep("diff.diff", nrule), rep("diff.same", nrule)), ntree)

rownames(sim.data) <-NULL

sim.data$topo <- as.factor(sim.data$topo)
sim.data$link.rule <- as.factor(sim.data$link.rule)

boxplot(sim.data$h2~sim.data$link.rule*sim.data$topo)
boxplot(sim.data$z.nodf[sim.data$z.nodf > 0]~sim.data$link.rule[sim.data$z.nodf > 0]*sim.data$topo[sim.data$z.nodf > 0])

nodf.aov <- aov(sim.data$z.nodf[sim.data$z.nodf > 0]~sim.data$link.rule[sim.data$z.nodf > 0]*sim.data$topo[sim.data$z.nodf > 0])
h2.aov <- aov(sim.data$h2~sim.data$link.rule*sim.data$topo)
ISA.aov <- aov(sim.data$z.ISA~sim.data$link.rule*sim.data$topo)

sim.data[which(sim.data$nofd)]

write.table(sim.data, file="sim_results_02_08.csv", sep=",")