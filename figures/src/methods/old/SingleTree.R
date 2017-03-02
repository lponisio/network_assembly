library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)
library(RColorBrewer)
setwd("~/Dropbox/network_assembly")
source('simulation/R/coal/SimTreeTrait_coal.R')
source('simulation/R/coal/initialize_coal.R')
prms <- base.prms()
prms$sp <- 20
prms$combinations <- expand.grid(plants=1:prms$sp,
                                 animals=1:prms$sp)
setseed(104)


plot.trees <- function(){
  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  path <- "~/Dropbox/network_assembly/figures/methods"
  all.trees <- sim.phylo(prms$sp)
  Traits <- character.evo(all.trees[[1]], all.trees[[2]])

  plot.trait <- function(tree, traits, dir, cols){
    standard.traits <- traits-min(traits)
    standard.traits <- standard.traits/max(standard.traits)
    plot(tree,show.tip.label=FALSE, direction=dir, edge.col=cols)
    tiplabels(pch=16,col=cols,cex=standard.traits*2)
  }

  f <- function(){
    layout(matrix(1:2, ncol=2))
    par(oma=c(2,2,2,2), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    plot.trait(tree=all.trees[[1]], traits=Traits$diff.same[[1]],
               dir="r", cols="white")
    plot.trait(tree=all.trees[[2]], traits=Traits$diff.same[[2]],
               dir="l", cols="white")
  }
  pdf.f(f, file= file.path(path, "trees.pdf"), width=5, height=5)
}

plot.trees()
