library(RColorBrewer)

diff.bar3 <- function(simres,
                      path.fig,
                      metric1,
                      metric2,
                      metric3,
                      subset=FALSE,
                      column=NA,
                      case=NA,
                      met1lab='Relative \n  nestness',
                      met2lab='Relative \n modularity',
                      met3lab='Connectance',
                      adj.lab=0.05){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }

  calc.diff.mean.sd <- function(met, link.type){
    sub.fun <- function(topo.type, met){
      simres[[1]][, met][simres[[2]][,'topo'] == topo.type  &
                         simres[[2]][,'link.rule'] == link.type]
    }
    coevo.difs <- function(mets){
      combins1 <- cbind(1:(length(mets[[1]])/2),
      ((length(mets[[1]])/2) + 1): length(mets[[1]]))
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
    prep.mets <-
      lapply(unique(simres[[2]][,'topo']), sub.fun, met=met)

    diff.mets <- coevo.difs(prep.mets)
    colnames(diff.mets) <- c('mean', 'sd', 'ci.lb', 'ci.ub')
    diff.mets <- diff.mets[!is.na(diff.mets[,'mean']),]
    return(diff.mets)
  }

  met.1.quan <- calc.diff.mean.sd(met=metric1,
                                  link.type='quan')
  met.1.qual <- calc.diff.mean.sd(met=metric1,
                                  link.type='qual')
  r.met1 <- range(met.1.quan, met.1.qual, na.rm=TRUE)
  met.2.quan <- calc.diff.mean.sd(met=metric2,
                                  link.type='quan')
  met.2.qual <- calc.diff.mean.sd(met=metric2,
                                  link.type='qual')
  r.met2 <- range(met.2.quan, met.2.qual, na.rm=TRUE)
  met.3.quan <- calc.diff.mean.sd(met=metric3,
                                  link.type='quan')
  met.3.qual <- calc.diff.mean.sd(met=metric3,
                                  link.type='qual')
  r.met3 <- range(met.3.quan, met.3.qual, na.rm=TRUE)

  out <- rbind(met.1.quan,  met.1.qual,
               met.2.quan,  met.2.qual,
               met.3.quan, met.3.qual)
  out <- apply(out, 2, round, 3)
  out <- data.frame(out)
  out$metric <- rep(c(metric1, metric2, metric3), each=4)
  out$weights <- rep(c("Unweighted", "Unweighted", "Weighted",
                       "Weighted"), 3)
  out$coevolution <- rep(c("with cospeciation",
                           "without cospeciation"), 6)
  out <- out[, c("coevolution", "metric",
                 "weights", "mean", "sd",
                 "ci.lb", "ci.ub")]
  out <- out[, c(-6, -7)]

  write.table(out, file=file.path(path.fig, "saved/diff.txt"),
              sep="&",
              row.names=FALSE)


  cols <-grey(c(0.5,0.8))
  f <- function(){
    plot.met <- function(sum.dats, r.qq, cols, labels=FALSE,...){
      mp <- barplot(sum.dats[,1], xlab='', ylab='',
                    col=cols,
                    las=1,
                    ylim=r.qq,
                    cex.axis=1.2,...)
      abline(h=0, col='gray')
      segments(mp, sum.dats[,2], mp, sum.dats[,3],
               lwd=2)
      if(labels){
        labels <- c('Effect of coevolution \n with cospeciation',
                    'Effect of coevolution \n without cospeciation')
        text(mp, par('usr')[3] - adj.lab,
             srt = 45, adj = 1,
             labels = labels,
             xpd = NA,
             cex=1)
      }
    }

    layout(matrix(1:6, 3, 2, byrow=FALSE),
           widths=rep(1,8), heights=rep(1,8))
    par(oma=c(7,3,2,0.1), mar=c(0.5,3,1.5,0.1),
        mgp=c(0,0.15,0), tcl=0, cex.axis=0.8)
    plot.met(met.1.qual, r.met1, cols, xaxt='n')
    mtext(met1lab, side=2, line=3, cex=1)
    mtext('Unweighted', side=3, line=1, cex=1.5)
    plot.met(met.2.qual, r.met2, cols, xaxt='n')
    mtext(met2lab, side=2, line=3, cex=1)
    plot.met(met.3.qual, r.met3, cols, labels=TRUE, xaxt='n')
    mtext(met3lab, side=2, line=3, cex=1)

    plot.met(met.1.quan, r.met1, cols, xaxt='n', yaxt='n')
    mtext('Weighted', side=3, line=1, cex=1.5)
    plot.met(met.2.quan, r.met2, cols, xaxt='n', yaxt='n')
    plot.met(met.3.quan, r.met3, cols, labels=TRUE,  xaxt='n',  yaxt='n')
  }

  pdf.f(f, file= file.path(path.fig,
                           sprintf('%s.pdf',
                                   paste('diff', metric1, metric2, metric3,
                                         sep=''))), height=6, width=5)

}
