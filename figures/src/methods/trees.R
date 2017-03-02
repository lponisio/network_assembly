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
source('simulation/src/bd/SimPhyloBD.R')
tree <- "bd"
source('simulation/src/all/initialize.R')
source('figures/src/misc.R')
prms <- base.prms()
prms$sp <- 20
prms$combinations <- expand.grid(plants=1:prms$sp,
                                 animals=1:prms$sp)


## shallow
cases.shallow <- c(mu= 0.01,
                   lambda= 0.05)

## deep
cases.deep <- c(mu= 0.5,
                lambda= 10)


set.seed(44)
plot.trees <- function(){
  deep <- sim.bd.taxa.age(n=prms$sp,
                          numbsim=1,
                          lambda=cases.deep["lambda"],
                          mu= cases.deep["mu"],
                          frac= 1,
                          age= prms$age,
                          mrca= FALSE)[[1]]

  deep.low <- rTraitCont(deep, model='BM', sigma=0.01)
  deep.high <- rTraitCont(deep, model='BM', sigma=4)

  shallow <- sim.bd.taxa.age(n=prms$sp,
                             numbsim=1,
                             lambda=cases.shallow["lambda"],
                             mu= cases.shallow["mu"],
                             frac= 1,
                             age= prms$age,
                             mrca= FALSE)[[1]]

  shallow.low <- rTraitCont(shallow, model='BM', sigma=0.01)
  shallow.high <- rTraitCont(shallow, model='BM', sigma=4)
  layout(matrix(c(1,1,2,3,4,4,
                  5,6),
                nrow=2 ,ncol=4),
         widths=c(3,3,3,3))
  par(oma=c(2, 2, 1, 1), mar=c(0.5, 1, 1, 1.5), mgp=c(2, 1, 0))
  plot(deep, show.tip.label=FALSE,
       main=expression(paste("Deep: ", mu,'=0.05, ', gamma, "=10")), cex=1.2)
  plot(density(deep.low), main="", yaxt="n", xaxt="n", type="l",
       lwd=2)
  legend("topleft", legend=expression(paste(sigma,'=0.01')), bty="n",
         cex=1.2)
  mtext("Frequency" , 2,line=1)

  plot(density(deep.high), main="", yaxt="n", xaxt="n", type="l", lwd=2)
  legend("topright", legend=expression(paste(sigma,'=4.0')), bty="n",
         cex=1.2)
  mtext("Trait values" , 1,line=1)
  mtext("Frequency" , 2,line=1)

  plot(shallow, show.tip.label=FALSE,
       main=expression(paste("Shallow: ", mu,'=0.01, ', gamma, "=0.05")), cex=1.2)
  plot(density(shallow.low), main="", yaxt="n", xaxt="n", type="l",
       lwd=2)
  mtext("Frequency" , 2,line=1)
  legend("topright", legend=expression(paste(sigma,'=0.01')), bty="n",
         cex=1.2)
  plot(density(shallow.high), main="", yaxt="n", xaxt="n", type="l",
       lwd=2)
  legend("topright", legend=expression(paste(sigma,'=4.0')), bty="n",
         cex=1.2)
  mtext("Trait values" , 1, line=1)
  mtext("Frequency" , 2,line=1)
}

path <- "figures/methods"
pdf.f(plot.trees, file=file.path(path, "trees.pdf"),
      height=5, width=7)

