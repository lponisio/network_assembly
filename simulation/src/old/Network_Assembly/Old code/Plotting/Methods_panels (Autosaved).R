library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)

setwd("~/Documents/Network_Assembly")
source("Functions/plot_trait_pch.R")
source("Functions/everything8.R")
#set.seed(22)
## first step
npol <- 15
ntree <- 1
all.trees <- sim.phylo(nsim=ntree, mu=0.1, lambda=0.1, age=1000, nspecies=npol)

Traits <- character.evo(all.trees[[1]][[1]], all.trees[[2]][[1]])

quartz(width=5,height=9)
layout(matrix(1:8, nrow=2, byrow=TRUE))
old.mar <- par("mar")

#layout(matrix(c(1,1,5,5,2,2,6,6,3,3,7,7,4,4,8,8,9,10,11,12), nrow=4,ncol=5),widths=c(rep(0.5,16), rep(1,4)))
par(mar = old.mar + c(-4,-2,-4, -2))

##same tree same traits
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.same[[1]], dir="r")
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="l")


## different tree, same traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$diff.same[[1]],  dir="r")
plot.trait(tree=all.trees[[2]][[1]], traits=Traits$diff.same[[2]],  dir="l")

## same tree, different traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="r")
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[2]],  dir="l")


## different tree, different traits

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$diff.diff[[1]],  dir="r")
plot.trait(tree=all.trees[[2]][[1]], traits=Traits$diff.diff[[2]],  dir="l")


## step 2

int.mats <- make.pa.mat(Traits$diff.diff)

bar <- int.mats$trait.bar1
bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
bar.cols <- ifelse(bar == 1, "darkred", 0)

comp <- int.mats$trait.narrow1
comp <- comp[order(rowSums(comp), decreasing=TRUE),]
comp <- comp[,order(colSums(comp), decreasing =TRUE)]

cInt.comp <- classIntervals(comp, style="fisher")
cPal <- tim.colors(24)
comp.cols  <- findColours(cInt.comp,cPal)
comp.cols <- ifelse(comp.cols == "#00008F", "white", comp.cols)

neu <- int.mats$neutral1

exp.int <- expand.grid(1:npol, 1:npol)

#layout(matrix(c(1,1,1, 2,2,2, 3,4,5), nrow=3, ncol=3), widths=c(0.5,0.5,1))
#layout.show(5)
quartz(width=3,height=9)
layout(matrix(1:3, nrow=3))
old.mar <- par("mar")
par(mar = old.mar + c(0,-1,0,-1))

#plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.same[[1]], dir="r")
#plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="l")

par(mar = old.mar +  c(-3,-1,-3,-1))
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(npol,1))
gradient.rect(xleft=14,xright=15, ytop=4, ybottom=1, col=unique(comp.cols[comp.cols!="white"]), gradient="y", border=NA)

plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(npol,1))
plot(exp.int[,1], exp.int[,2], col="darkred", pch=15, ylab="", xlab="", xaxt="n", yaxt="n")


#####################################################################
# with abund matrices too

abund.pol <- rlnorm(npol, meanlog=1, sdlog=1)
abund.plant <- rlnorm(npol, meanlog=1, sdlog=1)
abund <- outer(abund.plant, abund.pol)

trait.abund.comp <- comp*abund
trait.abund.bar <- bar*abund
trait.abund.neu <- neu*abund

abund.comp.cols <- prep.cols(trait.abund.comp)
abund.bar.cols <- prep.cols(trait.abund.bar)
abund.neu.cols <- prep.cols(trait.abund.neu)

layout(matrix(c(1,1,1, 2,2,2, 3,4,5,6,7,8), nrow=3, ncol=4), widths=c(0.5,0.5,1,1))
layout.show(8)

old.mar <- par("mar")
par(mar = old.mar + c(0,-1,0,-1))

plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.same[[1]], dir="r")
plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="l")

par(mar = old.mar +  c(-3,-1,-3,-1))
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,1], exp.int[,2], col="darkred", pch=15, ylab="", xlab="", xaxt="n", yaxt="n")

plot(exp.int[,1], exp.int[,2], col=abund.comp.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,2], exp.int[,1], col=abund.bar.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,1], exp.int[,2], col=abund.neu.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n")

#### with abundance being "multipled" to comp fundamental

cInt.abund <- classIntervals(abund, style="fisher")
abund.cols  <- findColours(cInt.abund,cPal)

cInt.abund.trait <- classIntervals(trait.abund, style="fisher")
abund.trait.cols  <- findColours(cInt.abund.trait,cPal)
abund.trait.cols  <- ifelse(abund.trait.cols  == "#00008F", "white", abund.trait.cols) 

layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2))
old.mar <- par("mar")
par(mar = old.mar + c(-3,-1,-3,-1))
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n")
plot(1:10,1:10, col="white", ylab="", xlab="", xaxt="n", yaxt="n", bty="n")
plot(1:10,1:10, col="white", ylab="", xlab="", xaxt="n", yaxt="n", bty="n")

plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n")
plot(exp.int[,1], exp.int[,2], col=abund.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n")
plot(exp.int[,1], exp.int[,2], col=abund.trait.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n")
## step 4

plotweb(trait.abund,
			method="normal",
			arrow="both",
			bor.col.interaction="grey",
			y.width.low=0.02,
			y.width.high=0.02, 
			col.high="darkolivegreen", 
			col.low="darkblue", 
			high.lablength=0, 
			low.lablength=0)


plotweb(comp,
			method="normal",
			arrow="both",
			bor.col.interaction="grey",
			y.width.low=0.02,
			y.width.high=0.02, 
			col.high="darkolivegreen", 
			col.low="darkblue", 
			high.lablength=0, 
			low.lablength=0)

			
rownames(comp) <- colnames(comp) <- 1:30
mod.comp <- computeModules(empty(round(comp*10)), steps=10^4)
mod.comp.web <- prepareWebForPlottingModules(mod.comp)

cInt.mod <- classIntervals(as.vector(mod.comp.web), style="fisher")
mod.cols  <- findColours(cInt.mod,cPal)
mod.cols <- ifelse(mod.cols  == "#00008F", "white", mod.cols) 
plot(exp.int[,1], exp.int[,2], col=mod.cols, pch=15, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
