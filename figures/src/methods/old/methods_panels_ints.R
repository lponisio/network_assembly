library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)
library(RColorBrewer)

setwd("~/Dropbox/network_assembly")
source('simulation/src/all/SimTreeTrait.R')
source('simulation/src/coal/all/SimPhyloCoal.R')
source('simulation/src/coal/varIntimacy_master.R')
source('simulation/src/coal/all/initialize.R')
source('figures/src/misc.R')
prms <- base.prms()
prms$sp <- 20
prms$combinations <- expand.grid(plants=1:prms$sp,
                                 animals=1:prms$sp)
set.seed(4444)


plot.methods <- function(){  
  path <- "figures/methods"
  all.trees <- sim.phylo(prms$sp)
  Traits <- character.evo(all.trees[[1]], all.trees[[2]])

  plot.trait <- function(tree, traits, dir, cols){
    standard.traits <- traits-min(traits)
    standard.traits <- standard.traits/max(standard.traits)
    plot(tree,show.tip.label=FALSE, direction=dir, edge.col=cols)
    tiplabels(pch=16, col=cols, cex=standard.traits*2)
  }

  ## prep.cols <- function(dats){
  ##   cInt <- classIntervals(dats, style="fisher")
  ##   cPal <-gray(100:0/100)
  ##   cols  <- findColours(cInt,cPal)
  ##   return(cols)
  ## }

  plot.comp <- function(int.mats, exp.int){
    comp <- int.mats$matching
    comp <- comp[order(rowSums(comp), decreasing=TRUE),]
    comp <- comp[,order(colSums(comp), decreasing =TRUE)]
    ## comp.cols <- prep.cols(comp)
    comp.cols <- ifelse(comp !=0, "black", 0)
    par(mar=c(0.5,2,5.5,2))
    plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15,
         ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1),
         cex=1.4)
    mtext("Matching", side=3, line=0.3, cex=1)
  }
  plot.barrier <- function(int.mats, exp.int){
    bar <- int.mats$barrior
    bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
    bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
    bar.cols <- ifelse(bar !=0, "black", 0)
    par(mar=c(2,2,4,2))
    plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
         xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1), cex=1.4)
    mtext("Barrier", side=3, line=0.3, cex=1)
  }
  
  f <- function(){
    layout(matrix(c(1,2,9,3,4,11,
                    1,2,10,3,4,12,
                    5,6,13,7,8,15,
                    5,6,14,7,8,16),
                  nrow=4 ,ncol=6, byrow=TRUE),
           widths=c(1,1,2,1,1,2))
    old.mar <- par("mar")
    par(mar = old.mar + c(-4,-1,-1, -2))
    topo.col <- brewer.pal(11, 'Spectral')[c(1,3,10,11)]
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    ##same tree same traits
    plot.trait(tree=all.trees[[1]], traits=Traits$same.same[[1]],
               dir="r", cols=topo.col["same.same"])
    mtext("Coevolution", side=2, line=0.5, cex=1.5)
    mtext("a)", side=3, adj=0, line=-0.5)
    mtext("Cospeciation", side=3, adj=-0.7, cex=1.5)
    plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[1]],
               dir="l", cols=topo.col["same.same"])
    ## different tree, same traits
    plot.trait(tree=all.trees[[1]], traits=Traits$diff.same[[1]],
               dir="r", cols=topo.col["diff.same"])
    mtext("No Cospeciation", side=3, adj=-0.4,  cex=1.5)
    mtext("b)", side=3, adj=0, line=-0.5)
    plot.trait(tree=all.trees[[2]], traits=Traits$diff.same[[2]],
               dir="l", cols=topo.col["diff.same"])
    ## same tree, different traits
    plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[1]],
               dir="r", cols=topo.col["same.diff"])
    mtext("No Coevolution", side=2, line=0.5,  cex=1.5)
    mtext("c)", side=3, adj=0, line=-0.5)
    plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[2]],
               dir="l", cols=topo.col["same.diff"])
    ## different tree, different traits
    plot.trait(tree=all.trees[[1]], traits=Traits$diff.diff[[1]],
               dir="r",cols=topo.col["diff.diff"])
    mtext("d)", side=3, adj=0, line=-0.5)
    plot.trait(tree=all.trees[[2]], traits=Traits$diff.diff[[2]],
               dir="l", cols=topo.col["diff.diff"])
    
    ##plot interaction matrices
    ss.mats <- make.pa.mat(prms, Traits$same.same[1:2], 1, 0.5)
    sd.mats <- make.pa.mat(prms, Traits$same.diff[1:2], 1, 0.5)
    ds.mats <- make.pa.mat(prms, Traits$diff.same[1:2], 1, 0.5)
    dd.mats <- make.pa.mat(prms, Traits$diff.diff[1:2], 1, 0.5)
    exp.int <- prms$combinations
    plot.comp(ss.mats, exp.int)
    plot.barrier(ss.mats, exp.int)
    plot.comp(sd.mats, exp.int)
    plot.barrier(sd.mats, exp.int)
    plot.comp(ds.mats, exp.int)
    plot.barrier(ds.mats, exp.int)
    plot.comp(dd.mats, exp.int)
    plot.barrier(dd.mats, exp.int)


    
  }
  pdf.f(f, file= file.path(path, "methods.pdf"), width=8, height=9)
}

plot.methods()

