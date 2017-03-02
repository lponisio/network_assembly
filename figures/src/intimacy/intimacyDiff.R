intInt.diff.plot3 <- function(simres,
                              metric1,
                              metric2,
                              metric3,
                              met1lab='Nestedness',
                              met2lab='Modularity',
                              met3lab,
                              xmetric,
                              xlabel,
                              leg.loc='bottomleft',
                              left.lab='Low',
                              right.lab='High',
                              diff.load=FALSE,
                              range.intint.quan,
                              range.intint.qual,
                              path.fig){

  calc.diff.mean.sd <- function(met, link.type, xmet){
    sub.fun <- function(topo.type, met, xmet){
      simres[[1]][, met][simres[[2]][,'topo'] == topo.type  &
                         simres[[2]][,'link.rule'] == link.type &
                         simres[[1]][, xmetric] == xmet]
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

    prep.mets <-lapply(
      unique(simres[[1]][, xmetric])[!is.na(unique(simres[[1]][, xmetric]))],
      function(x) {
        lapply(unique(simres[[2]][,'topo']), sub.fun, met=met, xmet=x)
      })
    min.reps <- rapply(prep.mets, function(x) sum(!is.na(x)),
                       how="list")
    to.keep <- sapply(min.reps, function(x) all(x > 10))
    prep.mets <- prep.mets[to.keep]
    min.reps <- min(rapply(prep.mets, function(x) sum(!is.na(x))))
    min.reps <- 2*round((min.reps/2) -0.5)
    prep.mets <- rapply(prep.mets, function(x) {
      x <- x[!is.na(x)]
      r.samp <- try(sample(1:length(x), min.reps))
      new.x <- x[r.samp]
      return(new.x)
    }, how="replace")

    diff.mets <- do.call(rbind, lapply(prep.mets, coevo.difs))
    diff.mets <- cbind(diff.mets,
                       rep(unique(
                         simres[[1]][,xmetric])[!is.na(unique(
                           simres[[1]][,  xmetric]))][to.keep],
                           each=2))
    colnames(diff.mets) <- c('mean', 'sd', 'ci.lb', 'ci.ub', 'intint')
    diff.mets <- diff.mets[order(diff.mets[,'intint']),]
    diff.mets <- diff.mets[!is.na(diff.mets[,'mean']),]
    return(diff.mets)
  }
  if(!diff.load){
    met.1.quan <- calc.diff.mean.sd(met=metric1,
                                    link.type='quan')
    met.1.qual <- calc.diff.mean.sd(met=metric1,
                                    link.type='qual')
    r.met1 <- range(met.1.quan, met.1.qual, na.rm=TRUE)  + c(-2, 2)
    met.2.quan <- calc.diff.mean.sd(met=metric2,
                                    link.type='quan')
    met.2.qual <- calc.diff.mean.sd(met=metric2,
                                    link.type='qual')
    r.met2 <- range(met.2.quan, met.2.qual, na.rm=TRUE) + c(0, 4)
    met.3.quan <- calc.diff.mean.sd(met=metric3,
                                    link.type='quan')
    met.3.qual <- calc.diff.mean.sd(met=metric3,
                                    link.type='qual')
    r.met3 <- range(met.3.quan, met.3.qual, na.rm=TRUE)
    save(met.1.quan,  met.1.qual,  met.2.quan,  met.2.qual,
         met.3.quan,  met.3.qual, r.met1, r.met2, r.met3,
         file='figures/src/intimacy/saved/diff.Rdata')
  } else  load(file='figures/src/intimacy/saved/diff.Rdata')

  ## cols <- brewer.pal(11, 'Spectral')[c(2,5)]
  cols <-grey(c(0.5,0.8))
  fill.col <- add.alpha(cols, alpha=0.2)

  f <- function(){
    plot.met <- function(dats, r.qq, r.intint,
                         cols,
                         fill.col,
                         rule,
                         plot.leg=FALSE,
                         lab=NA,...){
      plot(smooth.spline(dats[rownames(dats) == 'coevo.cosp','mean'] ~
                         dats[rownames(dats) == 'coevo.cosp', 'intint'],
                         df=5),
           xlab='', ylab='',
           col="black",
           las=1,
           ylim=r.qq,
           xlim=r.intint,
           yaxt='n',
           xaxt='n',
           cex.axis=1.2,
           type='l',
           lwd=2,...)

      points(smooth.spline(dats[rownames(dats) == 'coevo.cosp','ci.lb'] ~
                           dats[rownames(dats) == 'coevo.cosp', 'intint'],
                           df=5),
             col=cols[1],
             type='l',
             lty="solid",
             lwd=2)
      points(smooth.spline(dats[rownames(dats) == 'coevo.cosp','ci.ub'] ~
                           dats[rownames(dats) == 'coevo.cosp', 'intint'],
                           df=5),
             col=cols[1],
             type='l',
             lty="solid",
             lwd=2)

      points(smooth.spline(dats[rownames(dats) == 'coevo','mean'] ~
                           dats[rownames(dats) == 'coevo', 'intint'],
                           df=5),
             col="black",
             type='l',
             lty="dotted",
             lwd=2)

      points(smooth.spline(dats[rownames(dats) == 'coevo','ci.lb'] ~
                           dats[rownames(dats) == 'coevo', 'intint'],
                           df=5),
             col=cols[1],
             type='l',
             lty="dotted",
             lwd=2)
      points(smooth.spline(dats[rownames(dats) == 'coevo','ci.ub'] ~
                           dats[rownames(dats) == 'coevo', 'intint'],
                           df=5),
             col=cols[1],
              type='l',
             lty="dotted",
             lwd=2)


      out1.lb <- smooth.spline(
        dats[rownames(dats) == 'coevo.cosp', 'intint'],
        dats[rownames(dats) == 'coevo.cosp', 'ci.lb'], df=5)
      out1.ub <- smooth.spline(
        dats[rownames(dats) == 'coevo.cosp', 'intint'],
        dats[rownames(dats) == 'coevo.cosp', 'ci.ub'], df=5)
      out2.lb <- smooth.spline(
        dats[rownames(dats) == 'coevo', 'intint'],
        dats[rownames(dats) == 'coevo', 'ci.lb'], df=5)
      out2.ub <- smooth.spline(
        dats[rownames(dats) == 'coevo', 'intint'],
        dats[rownames(dats) == 'coevo', 'ci.ub'], df=5)
      y1.min <- out1.lb$y
      y1.max <-  out1.ub$y
      y2.min <- out2.lb$y
      y2.max <-  out2.ub$y

      polygon(x=c(dats[rownames(dats) == 'coevo.cosp', 'intint'],
                rev(dats[rownames(dats) == 'coevo.cosp', 'intint'])),
              y=c(y1.min, rev(y1.max)),
              col=fill.col[1], border = NA)
      polygon(x=c(dats[rownames(dats) == 'coevo', 'intint'],
                rev(dats[rownames(dats) == 'coevo', 'intint'])),
              y=c(y2.min, rev(y2.max)),
              col=fill.col[2], border = NA)
      ## abline(h=0, lty='solid', col='grey', lwd=2)
      if(rule == 'qual'){
        axis(2, at=pretty(r.qq, n=5, min.n=3),
             cex.axis=1.5, las=2)
        mtext(lab, side=2, line=5, cex=1.5)
      }
      if(plot.leg){
        legend(leg.loc,
               legend=c('Effect of coevolution with cospeciation',
                 'Effect of coevolution without cospeciation'),
               col="black", lty=c('solid', "dotted"), lwd=1.5, bty='n', cex=1.2)
      }
    }

    layout(matrix(1:6, ncol=2, byrow=TRUE))
    par(oma=c(8, 9, 5, 1),
        mar=c(0.5, 0.5, 0.5, 1),
        mgp=c(2, 1, 0))
    plot.met(met.1.qual, r.met1,
             range.intint.qual,
             cols,
             fill.col,
             rule="qual",
             lab=met1lab,
             plot.leg=TRUE)
    legend("topright", "a)", bty="n", cex=1.2)
    mtext('Unweighted', side=3, line=1, cex=1.5)
    plot.met(met.1.quan, r.met1,
             range.intint.quan,
             cols,
             fill.col,
             rule="quan")
    legend("topright", "b)", bty="n", cex=1.2)
    mtext('Weighted', side=3, line=1, cex=1.5)
    plot.met(met.2.qual, r.met2,
             range.intint.qual,
             cols,
             fill.col,
             rule="qual",
             lab=met2lab)
    legend("topright", "c)", bty="n", cex=1.2)
    plot.met(met.2.quan, r.met2,
             range.intint.quan,
             cols,
             fill.col,
             rule="quan")
    legend("topright", "d)", bty="n", cex=1.2)
    plot.met(met.3.qual, r.met3 + c(0,0.2),
             range.intint.qual,
             cols,
             fill.col,
             rule="qual",
             lab=met3lab)
    legend("topright", "e)", bty="n", cex=1.2)
    axis(1, at=pretty(range.intint.qual, n=5, min.n=3),
         cex.axis=1.5)
    mtext(xlabel, side=1, line=6.5, cex=1.5)
    mtext(left.lab, side=1, line=4, cex=1.5, adj=0)
    mtext(right.lab, side=1, line=4, cex=1.5, adj=1)
    plot.met(met.3.quan, r.met3 + c(0,0.2),
             range.intint.quan,
             cols,
             fill.col,
             rule="quan",
             plot.leg=FALSE)
    legend("topright", "f)", bty="n", cex=1.2)
    axis(1, at=pretty(range.intint.quan, n=5, min.n=3),
         cex.axis=1.5)
    mtext(xlabel, side=1, line=6.5, cex=1.5)
    mtext(left.lab, side=1, line=4, cex=1.5, adj=0)
    mtext(right.lab, side=1, line=4, cex=1.5, adj=1)
  }
  pdf.f(f, file= file.path(path.fig, sprintf('%s.pdf',
             paste('diff', xmetric,
                   metric1, metric2, metric3,
                   sep=''))), width=8, height=9)

}
