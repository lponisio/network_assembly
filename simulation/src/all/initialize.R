## ## packages to install
## install.packages("TreeSim")
## install.packages("ape")
## install.packages("geiger")
## install.packages("mvtnorm")
## install.packages("bipartite")
## install.packages("igraph")
## install.packages("parallel")
## install.packages("vegan")
## install.packages("RColorBrewer")

library(TreeSim)
library(ape)
library(geiger)
library(mvtnorm)
library(bipartite)
library(igraph)
library(parallel)
library(vegan)

## load source files
source('simulation/src/all/CalcMetrics.R')
source('simulation/src/all/LoadFix.R')
source('simulation/src/all/SimTreeTrait.R')

if(tree == 'coal'){
    ## use a coalescent process to generate trees
    source('simulation/src/coal/all/SimPhyloCoal.R')
    source('simulation/src/coal/CoalMaster.R')
    sim.dir <- 'simulation/saved/coal'
    fig.dir <- 'figures/coal'
    base.prms <- function() {
        prms.f(list(ntree= 1, ## do not change
                    nrule= 2,
                    nmat= 1,
                    nulls= 99,
                    sp= 15,
                    range.size=0.5,
                    sigma=0.1))
    }
} else if(tree == 'bd'){
    ## use a birth death process to generate trees
    source('simulation/src/bd/SimPhyloBD.R')
    source('simulation/src/bd/bdMaster.R')
    sim.dir <- 'simulation/saved/bd'
    fig.dir <- 'figures/bd'
    base.prms <- function() {
        prms.f(list(ntree= 1, ## do not change
                    nrule= 2,
                    nmat= 1,
                    nulls= 99,
                    sp= 15,
                    mu= 0.5,
                    lambda= 0.5,
                    age= 1000,
                    range.size=0.5,
                    sigma=0.1))
    }
}




prms.f <- function(prms) {
    prms$combinations <- expand.grid(plants=1:prms$sp,
                                     animals=1:prms$sp)
    ## everything that is created once and used throughout should be in
    ## here and then become a part of prms
    prms
}


## plotting code. Not included here because most was not used in the
## manuscript, can email me if you want it. lponisio@gmail.com

## source('figures/src/mod_by_nodf/mod_by_nodf_bar3.R')
## source('figures/src/mod_by_nodf/diff_bar3.R')
## source('figures/src/intimacy/intimacyDiff.R')
## source('figures/src/intimacy/intimacyDiff_image.R')
## source('figures/src/intimacy/intimacy2.R')
## source('figures/src/intimacy/intimacy3.R')
## source('figures/src/misc.R')
