library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
image.plot <- function(simres, metric, case){
  f <- function(){
    means.sim <- aggregate(list(metrics= simres[[1]][, metric]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"],
                             mu= unlist(simres[[1]][,"mu"]),
                             lambda= unlist(simres[[1]][,"lambda"])), 
                           function(x) mean(x, na.rm=TRUE))

    cats <- paste((means.sim$topo), (means.sim$link.rule),
                  (means.sim$mat))

    split.means <- split(means.sim, cats)

    prep.mats <- lapply(split.means, FUN=function(X){
      test <- matrix(X[,"metrics"], nrow=10, byrow=TRUE)
      rownames(test) <-unique(X[,"lambda"])
      colnames(test) <- unique(X[,"mu"])
      return(test)
    })

    all.dats <- range(means.sim$metrics[is.finite(means.sim$metrics)])  

    if(case=='evo'){
      to.plot <- prep.mats[grep("evo", names(prep.mats))]
    }
    if(case=='eco'){
      to.plot <- prep.mats[grep("eco", names(prep.mats))]
    }

    if(case=='eco2'){
      to.plot <- prep.mats[grep("eco2", names(prep.mats))]
    }

    matching <- to.plot[grep("matching", names(to.plot))]
    barrior  <- to.plot[grep("barrior", names(to.plot))]
    neutral <- to.plot[grep("neutral", names(to.plot))]

    to.plotz.side <- function(rule){
      lapply(rule, function(dats){
                                        #dats[dats < 0] <- 0
        ##print(dats)
        image(z=t(dats), x=as.numeric(colnames(dats)),
              y=as.numeric(rownames(dats)),
              ylab="", xlab="", tcl=0,
              col= gray(100:0/100), zlim= all.dats)
      }
             )
    }
    to.plotz <- function(rule){
      lapply(rule, function(dats) {
                                        #dats[dats < 0] <- 0
        ##print(dats)
        image(z=t(dats), x=as.numeric(colnames(dats)),
              y=as.numeric(rownames(dats)),
              ylab="", xlab="", tcl=0,
              col= gray(100:0/100),
              zlim= all.dats, yaxt="n")
      }
             )
    }
    layout(matrix(1:12, nrow= 4))
    par(las=1,cex.axis=1,cex.lab=1)
    par(oma=c(5,7,5,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    p <- to.plotz.side(matching)
    mtext(expression(lambda), side=2, line=3, outer=TRUE, las=0)
    mtext("Neutral evolution", side=2, line=5, outer=TRUE, las=0,
          adj=1)
    mtext('Coevo, no co-spec', side=2, line=5, outer=TRUE,
          las=0, adj=0.67)
    mtext('No Coevo, co-spec', side=2, line=5, outer=TRUE,
          las=0, adj=0.30)
    mtext('Coevo, co-spec', side=2, line=5, outer=TRUE,
          las=0, adj=-0.02)
    mtext('Matching', side=3, line=34)
    p <- to.plotz(barrior)
    mtext(expression(mu), side=1, line=3)
    mtext('Barrier', side=3, line=34)
    p <- to.plotz(neutral)
    mtext('Neutral', side=3, line=34)
  }

  path <- "~/Dropbox/network_assembly/figures/mu_lambda"
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(case, metric,
             sep=""))))
}

