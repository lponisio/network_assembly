##rm(list=ls())
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
source('simulation/src/coal/varAbund/initialize.R')
prms <- base.prms()
prms$sp <- 10
prms$combinations <- expand.grid(plants=1:prms$sp,
                                 animals=1:prms$sp)
##setseed(104)
setseed(104)

plot.link.rule <- function(){
  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  path <- "~/Dropbox/network_assembly/figures/methods"
  all.trees <- sim.phylo(prms$sp)
  Traits <- character.evo(all.trees[[1]], all.trees[[2]])

  prep.cols <- function(dats){
    cInt <- classIntervals(dats, style="fisher")
    cPal <-gray(0:100/100)
    cols  <- findColours(cInt,cPal)
    return(cols)
  }
  
  f <- function(){

    ##plot interaction matrices
    int.mats <- make.pa.mat(prms, Traits$diff.diff[1:2], 1, 0.5, 1)

    comp <- int.mats$matching.evo
    comp <- comp[order(rowSums(comp), decreasing=TRUE),]
    comp <- comp[,order(colSums(comp), decreasing =TRUE)]
    comp.cols <- prep.cols(comp)

    bar <- int.mats$barrior.evo
    bar <- bar[order(rowSums(bar),  decreasing=TRUE),]
    bar <- bar[,order(colSums(bar),  decreasing=TRUE)]
    bar.cols <- ifelse(bar != 0, comp.cols[1], 0)
    neu <- int.mats$neutral.evo

    exp.int <- prms$combinations
    layout(matrix(1:3, ncol=3))
    par(oma=c(4,4,4,4), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    plot(exp.int[,1], exp.int[,2], col=comp.cols, pch=15,
         ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1),
         cex=1.8)
    mtext("Matching", side=3, line=0.3, col="white")
    plot(exp.int[,2], exp.int[,1], col=bar.cols, pch=15, ylab="",
         xlab="", xaxt="n", yaxt="n", ylim=c(prms$sp,1), cex=1.8)
    mtext("Barrier", side=3, line=0.3,  col="white")
    plot(exp.int[,1], exp.int[,2], col=comp.cols[1], pch=15, ylab="",
         xlab="", xaxt="n", yaxt="n", cex=1.8)
    mtext("Neutral", side=3, line=0.3, col="white")
  }
  pdf.f(f, file= file.path(path, "linkrule.pdf"), width=9, height=3)
}

plot.link.rule()

