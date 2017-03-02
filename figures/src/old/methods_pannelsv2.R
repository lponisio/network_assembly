library(mvtnorm)
library(classInt)
library(fields)
library(TreeSim)
library(ape)
library(geiger)
library(plotrix)
setwd("~/Dropbox/network_assembly")
source('simulation/R/coal/SimTreeTrait_coal.R')
source('simulation/R/coal/initialize_coal.R')
prms <- base.prms()
set.seed(100)


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

  plot.trees <- function(){
    plot.trait <- function(tree, traits, dir){
      standard.traits <- traits-min(traits)
      standard.traits <- standard.traits/max(standard.traits)
      plot(tree,show.tip.label=FALSE, direction=dir)
      tiplabels(pch=16,col="black",cex=standard.traits*3)
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
      #layout(matrix(1:8, nrow=2, byrow=TRUE))
      old.mar <- par("mar")
      par(mar = old.mar + c(-4,-2,-2, -2))
      ##same tree same traits
      plot.trait(tree=all.trees[[1]], traits=Traits$same.same[[1]],
                 dir="r")
      mtext("Coevolution", side=2)
      mtext("Co-speciation", side=3, adj=-8)
      plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[1]],
                 dir="l")
      ## different tree, same traits
      plot.trait(tree=all.trees[[1]], traits=Traits$diff.same[[1]],
                 dir="r")
      mtext("No Co-speciation", side=3, adj=-1)
      plot.trait(tree=all.trees[[2]], traits=Traits$diff.same[[2]],
                 dir="l")
      ## same tree, different traits
      plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[1]],
                 dir="r")
      mtext("No Coevolution", side=2)
      plot.trait(tree=all.trees[[1]], traits=Traits$same.diff[[2]],
                 dir="l")
      ## different tree, different traits
      plot.trait(tree=all.trees[[1]], traits=Traits$diff.diff[[1]],
                 dir="r")
      plot.trait(tree=all.trees[[2]], traits=Traits$diff.diff[[2]],
                 dir="l")
    }
    #pdf.f(f, file= file.path(path, "trees.pdf"), width=5, height=9)

  }

  plot.mats <- function(){
    
    prep.cols <- function(dats){
      cInt <- classIntervals(dats, style="fisher")
      cPal <- tim.colors(24)
      cols  <- findColours(cInt,cPal)
      cols   <- ifelse(cols == "#00008F", "white", cols ) 
      return(cols)
    }

    f <- function(){
      int.mats <- make.pa.mat(prms, Traits$diff.diff[1:2], 1, 0.5)

      comp <- int.mats$matching.evo
      comp <- comp[order(rowSums(comp), decreasing=TRUE),]
      comp <- comp[,order(colSums(comp), decreasing =TRUE)]

      cInt.comp <- classIntervals(comp, style="fisher")
      cPal <- tim.colors(24)
      comp.cols  <- findColours(cInt.comp,cPal)
      comp.cols <- ifelse(comp.cols == "#00008F", "white", comp.cols)

      bar <- int.mats$barrior.evo
      bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
      bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
      bar.cols <- ifelse(bar == 1, comp.cols[1], 0)

      neu <- int.mats$neutral.evo

      exp.int <- prms$combinations

      #layout(matrix(1:3, nrow=3))
      old.mar <- par("mar")
      par(mar = old.mar + c(0,-1,0,-1))

      par(mar = old.mar +  c(-4,-1,-2,-1))
      plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15,
           ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1),
           cex=1.5)
      mtext("Matching", side=3, line=0.3)

      ##gradient.rect(xleft=14, xright=15, ytop=2, ybottom=1,
      ##              col=unique(comp.cols[comp.cols!="white"]),
      ##              gradient="y", border=NA)

      plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
           xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1), cex=1.5)
      mtext("Barrior", side=3, line=0.3)
      plot(exp.int[,1], exp.int[,2], col=comp.cols[1], pch=15, ylab="",
           xlab="", xaxt="n", yaxt="n", cex=1.5)
      mtext("Neutral", side=3, line=0.3)
    }
    pdf.f(f, file= file.path(path, "mats.pdf"), width=8, height=9)

  }
  plot.trees()
  plot.mats()
}
