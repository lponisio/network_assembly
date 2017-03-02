library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)
setwd("~/Dropbox/network_assembly")
source("figures/code/plot_trait_pch.R")
source("simulation/R/coal/SimTreeTrait_coal.R")
set.seed(22)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

prep.cols <- function(dats){
  cInt <- classIntervals(dats, style="fisher")
  cPal <- tim.colors(24)
  cols  <- findColours(cInt,cPal)
  cols   <- ifelse(cols == "#00008F", "white", cols ) 
  return(cols)
}

plot.trait <- function(tree, traits, dir){
  standard.traits <- traits-min(traits)
  standard.traits <- standard.traits/max(standard.traits)
  plot(tree,show.tip.label=FALSE, direction=dir)
  tiplabels(pch=16,col="black",cex=standard.traits*3)
}

plot.methods <- function(){
  all.trees <- sim.phylo(n=15)

Traits <- character.evo(all.trees[[1]][[1]], all.trees[[2]][[1]])

quartz(width=5,height=9)
layout(matrix(1:8, nrow=2, byrow=TRUE))
old.mar <- par("mar")

##layout(matrix(c(1,1,5,5,2,2,6,6,3,3,7,7,4,4,8,8,9,10,11,12), nrow=4,ncol=5),widths=c(rep(0.5,16), rep(1,4)))
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

int.mats <- make.pa.mat(Traits$diff.diff, 1, 0.5)

bar <- int.mats$barrior.evo
bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
bar.cols <- ifelse(bar == 1, "darkred", 0)

comp <- int.mats$matching.evo
comp <- comp[order(rowSums(comp), decreasing=TRUE),]
comp <- comp[,order(colSums(comp), decreasing =TRUE)]

cInt.comp <- classIntervals(comp, style="fisher")
cPal <- tim.colors(24)
comp.cols  <- findColours(cInt.comp,cPal)
comp.cols <- ifelse(comp.cols == "#00008F", "white", comp.cols)

neu <- int.mats$neutral.evo

exp.int <- expand.grid(1:npol, 1:npol)

##layout(matrix(c(1,1,1, 2,2,2, 3,4,5), nrow=3, ncol=3), widths=c(0.5,0.5,1))
##layout.show(5)
quartz(width=3,height=9)
layout(matrix(1:3, nrow=3))
old.mar <- par("mar")
par(mar = old.mar + c(0,-1,0,-1))

##plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.same[[1]], dir="r")
##plot.trait(tree=all.trees[[1]][[1]], traits=Traits$same.diff[[1]],  dir="l")

par(mar = old.mar +  c(-3,-1,-3,-1))
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15,
     ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(npol,1), cex=1.5)

gradient.rect(xleft=14,xright=15, ytop=4, ybottom=1,
              col=unique(comp.cols[comp.cols!="white"]),
              gradient="y", border=NA)

plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", ylim=c(npol,1), cex=1.5)

plot(exp.int[,1], exp.int[,2], col="darkred", pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", cex=1.5)


#####################################################################
                                        # with abund matrices too

trait.abund.bar <- int.mats$barrior.eco
trait.abund.comp <- int.mats$matching.eco
trait.abund.neu <- int.mats$neutral.eco

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
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,1], exp.int[,2], col="darkred", pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")

plot(exp.int[,1], exp.int[,2], col=abund.comp.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,2], exp.int[,1], col=abund.bar.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n", ylim=c(30,1))
plot(exp.int[,1], exp.int[,2], col=abund.neu.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")

#### with abundance being "multipled" to comp fundamental

cInt.abund <- classIntervals(abund, style="fisher")
abund.cols  <- findColours(cInt.abund,cPal)

cInt.abund.trait <- classIntervals(trait.abund, style="fisher")
abund.trait.cols  <- findColours(cInt.abund.trait,cPal)
abund.trait.cols  <- ifelse(abund.trait.cols  == "#00008F", "white", abund.trait.cols) 

layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2))
old.mar <- par("mar")
par(mar = old.mar + c(-3,-1,-3,-1))
plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")
plot(1:10,1:10, col="white", ylab="",
     xlab="", xaxt="n", yaxt="n", bty="n")
plot(1:10,1:10, col="white", ylab="",
     xlab="", xaxt="n", yaxt="n", bty="n")

plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")
plot(exp.int[,1], exp.int[,2], col=abund.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")
plot(exp.int[,1], exp.int[,2], col=abund.trait.cols, pch=15, ylab="",
     xlab="", xaxt="n", yaxt="n")
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
