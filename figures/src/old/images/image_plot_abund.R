library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
image.plot <- function(simres, metric){
  prep.cols <- function(dats){
    cInt <- classIntervals(dats, style="fisher")
    cPal <- tim.colors(10)
    cols  <- findColours(cInt,cPal)
    ##cols   <- ifelse(cols == "#00008F", "white", cols ) 
    return(cols)
  }
  f <- function(){
    means.sim <- aggregate(list(metrics= simres[[1]][, metric]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"],
                             abund.sd= unlist(simres[[1]][,"sd.abund"]),
                             abund.mean=unlist(simres[[1]][,"mean.abund"])),
                           function(x) mean(x, na.rm=TRUE))
    cats <- paste((means.sim$topo), (means.sim$link.rule),
                  (means.sim$mat))
    split.means <- split(means.sim, cats)
    
    prep.mats <- lapply(split.means, FUN=function(X){
      out.mat <- matrix(X[,"metrics"],
                        nrow=length(unique(X[,"abund.sd"])), byrow=TRUE)
      rownames(out.mat) <- unique(X[,"abund.sd"])
      colnames(out.mat) <- round(unique(X[,"abund.mean"]),2)
      return(out.mat)
    }) 
    
    to.plot <- prep.mats[grep("eco", names(prep.mats))]
    all.dats <- range(unlist(to.plot), na.rm=TRUE)
    matching <- to.plot[grep("matching", names(to.plot))]
    barrior  <- to.plot[grep("barrior", names(to.plot))]
    neutral <- to.plot[grep("neutral", names(to.plot))]

    to.plotz.side <- function(rule){
      lapply(rule, function(dats){
        #dats[dats < 0] <- 0
        print(dats)
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
        print(dats)
        image(z=t(dats), x=as.numeric(colnames(dats)),
              y=as.numeric(rownames(dats)),
              ylab="", xlab="", tcl=0,
              col= gray(100:0/100),
              zlim= all.dats, yaxt="n")
      }
             )
    }
    ## quartz()
    layout(matrix(1:12, nrow= 4))
    par(las=1,cex.axis=1,cex.lab=1)
    par(oma=c(5,7,5,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    p <- to.plotz.side(matching)
    mtext(expression(sigma), side=2, line=3, outer=TRUE, las=0)
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
  path <- "~/Dropbox/network_assembly/figures/coal/abund"
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("abund", metric,
             sep=""))))

}

