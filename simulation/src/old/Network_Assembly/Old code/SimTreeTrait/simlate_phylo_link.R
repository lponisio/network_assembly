
source("/Users/laurenponisio/Documents/Network Assembly /master_fun.R")
source("/Users/laurenponisio/Documents/Network Assembly /linkage_rules_faster.R")
source("/Users/laurenponisio/Documents/Network Assembly /simulate_character_evo.R")
source("/Users/laurenponisio/Documents/Network Assembly /simulate_trees.R")
source("/Users/laurenponisio/Documents/Network Assembly /all_network_mets.R")

Rprof()
set.seed(1)

ntree <- 2
ntopo <- 4
nrule <- 5

all.trees <- sim.phylo(nsim=ntree, mu=0.5, lambda=0.5, age=1000, nspecies=20)

sim.links <- mapply(master.fun,tree1=all.trees[[1]],tree2=all.trees[[2]], SIMPLIFY=FALSE)

sim.stats <-rapply(sim.links, f=network.metrics, how="replace")

sim.data <- t(matrix(rapply(sim.stats, c),nrow=9))

sim.data<- cbind(rep(c("trait.comp", "trait.bar", "neutral", "trait.narrow", "trait.wide"), ntopo*ntree), sim.data)
	
	
sim.data <- cbind(rep(c(rep("same.same", nrule), rep("same.diff", nrule), rep("diff.diff", nrule), rep("diff.same", nrule)), ntree), sim.data)
	
	
colnames(sim.data) <- c( "topo", "linkr", "nofd", "pval.nodf", "ISA", "pval.ISA", "h2", "pval.h2", "mod", "pval.mod", "num.mod")
	
Rprof(NULL)
summaryRprof()