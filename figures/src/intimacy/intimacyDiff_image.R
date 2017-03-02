library(fields)
intInt.diff.image <- function(simres,
                              metric1,
                              xmetric1,
                              xmetric2,
                              path.fig){
  
  calc.diff.mean.sd <- function(met, link.type, xmet){
    sub.fun <- function(topo.type, met, xmet){
      simres[[1]][, met][simres[[2]][,'topo'] == topo.type  &
                         simres[[2]][,'link.rule'] == link.type &
                         simres[[1]][, xmetric1] == xmet[1] &
                         simres[[1]][, xmetric2] == xmet[2]]
    }
    coevo.difs <- function(mets){
      combins1 <- cbind(1:(length(mets[[1]])/2),
                        ((length(mets[[1]])/2) + 1):
                        length(mets[[1]]))
      ## effect of coevolution with cospeciation
      mets.coevo.cosp <- mets[[1]][combins1[,1]] - 
        mets[[2]][combins1[,2]]

      ## effect of coevolution, no cospeciation
      mets.coevo <- mets[[4]][combins1[,1]] - 
        mets[[3]][combins1[,2]]
      
      sum.mets.coevo.cosp <-  mean.sd(mets.coevo.cosp)
      sum.mets.coevo <-  mean.sd(mets.coevo)
      return(rbind(coevo.cosp=sum.mets.coevo.cosp,
                   coevo=sum.mets.coevo))
    }
    ## unique combinations of specialization and partner overlap
    x1 <-  sort(unique(simres[[1]][, xmetric1])[!is.na(unique(simres[[1]][,
                                                                          xmetric1]))])
    x2 <- sort(unique(simres[[1]][, xmetric2])[!is.na(unique(simres[[1]][,
                                                                         xmetric2]))])
    combins <- expand.grid(x1,  x2)
    prep.mets <-apply(
      combins, 1,
      function(x) {
        lapply(unique(simres[[2]][,'topo']), sub.fun, met=met, xmet=x)
      })
    min.reps <- rapply(prep.mets, function(x) sum(!is.na(x)),
                       how="list")
    to.keep <- sapply(min.reps, function(x) all(x > 2))
    prep.mets <- prep.mets[to.keep]
    min.reps <- min(rapply(prep.mets, function(x) sum(!is.na(x))))
    min.reps <- 2*round((min.reps/2) -0.5)
    prep.mets <- rapply(prep.mets, function(x) {
      x <- x[!is.na(x)]
      r.samp <- try(sample(1:length(x), min.reps))
      new.x <- x[r.samp]
      return(new.x)
    }, how="replace")

    diff.mets <-  lapply(prep.mets, coevo.difs)
    diff.coevo.cosp <- cbind(do.call(rbind,
                                     lapply(diff.mets, function(x) x[1,])),
                             combins[to.keep,])
    diff.coevo <- cbind(do.call(rbind,
                                lapply(diff.mets, function(x) x[2,])),
                        combins[to.keep,])

    colnames(diff.coevo.cosp) <- colnames(diff.coevo) <-
      c('mean', 'sd', 'ci.lb', 'ci.ub', xmetric1, xmetric2)
    make.image.mat <- function(diff.mat){
      image.mat <- matrix(NA, ncol=length(x1), nrow=length(x2))
      colnames(image.mat) <- x1
      rownames(image.mat) <- x2
      image.mat[cbind(match(diff.mat[, xmetric2],
                            rownames(image.mat)),
                      match(diff.mat[, xmetric1],
                            colnames(image.mat)))] <-
                              diff.mat$mean
      return(image.mat)
    }
    coevo.cosp.image <- make.image.mat(diff.coevo.cosp)
    coevo.image <- make.image.mat(diff.coevo)
    all.vals <- range(diff.coevo$mean, diff.coevo.cosp$mean,
                      na.rm=TRUE)
    
    return(list(coevo.cosp= coevo.cosp.image,
                coevo = coevo.image,
                all.vals=all.vals))
  }
  met.quan <- calc.diff.mean.sd(met=metric1,
                                link.type='quan')
  met.qual <- calc.diff.mean.sd(met=metric1,
                                link.type='qual')

  make.image <- function(dats, all.dats, FUN = image.plot, ...){
    FUN(z=t(dats), x=as.numeric(colnames(dats)),
        y=as.numeric(rownames(dats)),
        ylab="", xlab="", tcl=0,
        col= brewer.pal(11, "RdYlGn"),
        zlim= all.dats,
        las=1, ...)
  }
  plot.image <- function(){
    make.image <- function(dats, all.dats, FUN = image, ...){
      FUN(z=t(dats), x=as.numeric(colnames(dats)),
          y=as.numeric(rownames(dats)),
          ylab="", xlab="", tcl=0,
          col= brewer.pal(11, "RdYlGn"),
          zlim= all.dats,
          las=1, ...)
    }
    layout(matrix(1:4, ncol=2, byrow=TRUE))
    par(oma=c(7, 9, 5, 3),
        mar=c(0.5, 0.5, 0.5, 1),
        mgp=c(2, 1, 0), cex.axis=1.5)
    make.image(met.quan$coevo.cosp, range(met.qual$all.vals,
                                     met.quan$all.vals),
               xaxt="n")
    mtext("Partner overlap", 2, line=4.5, cex=1.5, adj=-4.5)
    mtext("Unweighted", 3, line=2, cex=1.5)
    mtext("Weighted", 3, line=2, adj=2.5, cex=1.5) 
    make.image(met.quan$coevo.cosp, range(met.quan$all.vals,
                                     met.quan$all.vals),
               xaxt="n", yaxt="n",
               FUN=image.plot,
               legend.width=2)
    make.image(met.quan$coevo, range(met.quan$all.vals,
                                          met.quan$all.vals))
    mtext("Specialization", 1, cex=1.5, line=4, adj=3)
    make.image(met.quan$coevo, range(met.quan$all.vals,
                                          met.quan$all.vals),
               yaxt="n")
  }

  pdf.f(plot.image, file= file.path(path.fig, sprintf('%s.pdf',
                      paste("image", metric1,
                            sep=''))), width=7, height=6)

}
