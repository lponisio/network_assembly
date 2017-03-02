library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)
library(RColorBrewer)
setwd("~/Dropbox/network_assembly")
source('simulation/src/all/SimTreeTrait_samp_mods.R')
source('simulation/src/coal/SimPhyloCoal.R')
source('simulation/src/coal/varAbund/master.R')
source('simulation/src/coal/varAbund/initialize.R')
prms <- base.prms()
prms$sp <- 20
prms$combinations <- expand.grid(plants=1:prms$sp,
                                   animals=1:prms$sp)
set.seed(1)


plot.methods <- function(){
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

  prep.cols <- function(dats){
    cInt <- classIntervals(dats, style="fisher")
    cPal <-gray(100:0/100)
    cols  <- findColours(cInt,cPal)
    return(cols)
  }
  
  f <- function(){
    
    layout(matrix(c(1,2,3,4,9,
                    1,2,3,4,9,
                    1,2,3,4,10,
                    5,6,7,8,10,
                    5,6,7,8,11,
                    5,6,7,8,11),
                  nrow=6,ncol=5, byrow=TRUE),
           widths=c(rep(1,4), 2))

    old.mar <- par("mar")
    par(mar = old.mar + c(-4,-1,-1, -2))
    topo.col <- brewer.pal(4, 'Spectral')
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
    int.mats <- make.pa.mat(prms, Traits$diff.diff[1:2], 1, 0.5, 1)

    comp <- int.mats$matching.evo
    comp <- comp[order(rowSums(comp), decreasing=TRUE),]
    comp <- comp[,order(colSums(comp), decreasing =TRUE)]

    comp.cols <- prep.cols(comp)

    bar <- int.mats$barrior.evo
    bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
    bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
    bar.cols <- ifelse(bar !=0, comp.cols[1], 0)

    neu <- int.mats$neutral.evo

    exp.int <- prms$combinations
    par(mar=c(3,4,3,2))
    plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15,
         ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1),
         cex=1.8)
    mtext("Matching", side=3, line=0.3, cex=1.5)
    mtext("e)", side=3, line=0.5, adj=0)

    ##gradient.rect(xleft=14, xright=15, ytop=2, ybottom=1,
    ##              col=unique(comp.cols[comp.cols!="white"]),
    ##              gradient="y", border=NA)
    par(mar=c(3,4,3,2))
    plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
         xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1), cex=1.8)
    mtext("Barrier", side=3, line=0.3, cex=1.5)
    mtext("f)", side=3, line=0.5, adj=0)
    plot(exp.int[,1], exp.int[,2], col=comp.cols[1], pch=15, ylab="",
         xlab="", xaxt="n", yaxt="n", cex=1.8)
    mtext("Neutral", side=3, line=0.3, cex=1.5)
    mtext("g)", side=3, line=0.5, adj=0)
  }
  pdf.f(f, file= file.path(path, "methods.pdf"), width=8, height=9)
}

plot.methods()

