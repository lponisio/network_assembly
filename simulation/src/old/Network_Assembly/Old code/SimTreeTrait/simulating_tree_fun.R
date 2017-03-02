rm(list=ls())
library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)

x <-sample(1:1000,1)
x
set.seed(x)
tree1 <- sim.bd.taxa.age(n=15, numbsim=1, lambda=0.1, mu=0.1, frac=1, age=1000, mrca=TRUE)
tree2 <- sim.bd.taxa.age(n=15, numbsim=1, lambda=0.5,
                            mu=0.5, frac=1, age=1000, mrca=TRUE)

tree1.trait1 <- exp(rTraitCont(tree1[[1]], model="BM"))
sample.trait <- sample(log(tree1.trait1), size=Ntip(tree2[[1]]), replace=TRUE)
tree2.trait2 <-
    as.vector(exp(rmvnorm(1, sample.trait, sigma=0.01*vcv(tree2[[1]]),
                          method="chol")))
plot(tree1.trait1, tree2.trait2);abline(0,1)
                                    
tree1.trait1 <- rTraitCont(tree1[[1]], model="BM")
sample.trait <- sample(tree1.trait1, size=Ntip(tree2[[1]]), replace=TRUE)
tree2.trait2 <-
    as.vector(rmvnorm(1, sample.trait, sigma=0.01*vcv(tree2[[1]]),
                          method="chol"))

plot(exp(tree1.trait1), exp(tree2.trait2));abline(0,1)


##

tree1.trait1 <- exp(rTraitCont(tree1[[1]], model="BM"))

sample.trait <- sample(log(tree1.trait1), size=Ntip(tree2[[1]]), replace=TRUE)

names(sample.trait) <- tree2[[1]]$tip.label

tree2[[1]] <- reorder(tree2[[1]], "p")

alltraits.tree2 <- c(sample.trait,ace(x=sample.trait, phy=tree2[[1]])$ace)
alltraits.tree2 <- alltraits.tree2[tree2[[1]]$edge[,2]]

tree2.trait2 <- exp(rTraitCont(tree2[[1]], model="OU", theta=alltraits.tree2))

plot(density(tree1.trait1))
lines(density(tree2.trait2))


bla.sim <- rmvnorm(1,bla,sigma=0.05*vcv(tre),method="chol")
plot(bla,bla.sim);abline(0,1)


plot(tre)

bla.ace <- ace(bla.sim,tre)
cInt <- classIntervals(bla.ace$ace,24,style="fisher")
cPal <- tim.colors(24)
plot(cInt,cPal)
these.col <- findColours(cInt,cPal)

plot(tre,show.tip.label=FALSE)
nodelabels(pch=16,col=these.col,cex=1.5)


bla.sim2 <- rmvnorm(1,sigma=vcv(tre))
bla.ace2 <- ace(bla.sim2,tre)

cInt <- classIntervals(bla.ace2$ace,24,style="fisher")
cPal <- tim.colors(24)
plot(cInt,cPal)
these.col <- findColours(cInt,cPal)

plot(tre,show.tip.label=FALSE)
nodelabels(pch=16,col=these.col,cex=1.5)


####

all.trees <- sim.phylo(nsim=2, mu=0.5, lambda=0.5, age=1000, nspecies=15)

Traits <- character.evo(all.trees[[1]][[1]], all.trees[[2]][[1]])


##same tree same traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.same[[1]], dir="r")
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="l")

## same tree, different traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="r")
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[2]],  dir="l")


## different tree, different traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$diff.diff[[1]],  dir="r")
plot.trait(tree=all.trees[[2]][[1]], traits=Traits$diff.diff[[2]],  dir="l")


## different tree, same traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$diff.same[[1]],  dir="r")
plot.trait(tree=all.trees[[2]][[1]], traits=Traits$diff.same[[2]],  dir="l")

plot.trait <- function(tree, traits, dir){
	cInt <- classIntervals(traits,24,style="fisher")
	cPal <- tim.colors(24)
	#plot(cInt,cPal)
	these.col <- findColours(cInt,cPal)

	plot(tree,show.tip.label=FALSE, direction=dir)
	tiplabels(pch=16,col=these.col,cex=1.5)

}