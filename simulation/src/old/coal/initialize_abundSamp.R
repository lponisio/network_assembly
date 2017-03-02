library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)
library(parallel)
library(poilog)
library(moments)

base.prms <- function() {
  prms.f(list(ntree= 1, ##do not change
              nrule= 3,
              nmat= 11,
              nulls= 99,
              sp= 30,
              sd.abund= 1,
              mean.abund= 1,
              range.size=0.5))
}

prms.f <- function(prms) {
  prms$combinations <- expand.grid(plants=1:prms$sp,
                                   animals=1:prms$sp)
  ## everything that is created once and used throughout should be in
  ## here and then become a part of prms
  prms
}


