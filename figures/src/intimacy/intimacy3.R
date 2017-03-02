intInt.plot3 <- function(simres,
                         metric1,
                         metric2,
                         metric3,
                         met1lab='Nestedness',
                         met2lab='Modularity',
                         met3lab,
                         xmetric,
                         xlabel,
                         path.dir,
                         subset=FALSE,
                         column=NA,
                         case=NA,
                         leg.panel='right',
                         leg.loc='bottomleft',
                         rev.x,
                         range.intint=NA,
                         left.lab='Low',
                         right.lab='High',
                         leg.row=1){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }

  ## metric 1
  calc.stats <- function(x){
    m <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)
    ci.ub <- m + s
    ci.lb <- m - s
    return(c(mean=m, sd=s, ci.ub=ci.ub, ci.lb=ci.lb))
  }

  res.sim <- aggregate(list(met1= simres[[1]][, metric1],
                            met2= simres[[1]][, metric2],
                            met3= simres[[1]][, metric3]),
                       by= list(topo= simres[[2]][,'topo'],
                         link.rule= simres[[2]][,'link.rule'],
                         mat= simres[[2]][,'mats'],
                         intint = simres[[1]][, xmetric]),
                       calc.stats)
  ## metric 1
  range.met1 <- range(c(res.sim$met1[,'ci.lb'],
                        res.sim$met1[,'ci.ub']),
                      na.rm=TRUE)
  ## metric 2
  range.met2 <- range(c(res.sim$met2[,'ci.lb'],
                        res.sim$met2[,'ci.ub']),
                      na.rm=TRUE)
  ## metric 3
  range.met3 <- range(c(res.sim$met3[,'ci.lb'],
                        res.sim$met3[,'ci.ub']),
                      na.rm=TRUE)

  res.evo <- res.sim[res.sim$mat == 'evo',]
  res.evo <- res.evo[!apply(res.evo, 1, function(x) any(is.na(x))),]

  if(is.na(range.intint[1])) {
    if(rev.x) range.intint <- rev(range(res.evo$intint, na.rm =TRUE))
    else  range.intint <- range(res.evo$intint, na.rm =TRUE)
  }

  topo.col <- brewer.pal(11, 'Spectral')[c(1,3,10,11)]
  fill.col <- add.alpha(topo.col, alpha=0.2)
  names(topo.col) <- names(fill.col) <-
    c('diff.diff','same.diff','diff.same','same.same')
  f <- function(){
    plot.met1 <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        if(topos == 'same.same'){
          plot(smooth.spline(met1[,'mean'] ~ intint, df=5), xlab='', ylab='',
               col=topo.col[topo],
               las=1,
               ylim=range.met1,
               xlim=range.intint,
               yaxt='n',
               xaxt='n',
               cex.axis=1.2,
               type='l',
               lwd=2,...)
        }
        not.na <- !is.na(met1[,'ci.ub'])
        out1 <- smooth.spline(intint[not.na], met1[,'ci.lb'][not.na], df=5)
        out2 <- smooth.spline(intint[not.na], met1[,'ci.ub'][not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y
        c(y1.min, rev(y1.max))

        polygon(x=c(intint[not.na], rev(intint[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)

        points(smooth.spline(met1[,'mean'] ~ intint, df=5),
               col=topo.col[topo],
               type='l', lwd=2, ...)
        if(topos == 'same.same'){
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(rule == 'qual'){
            axis(2, at=pretty(range.met1, n=5, min.n=3),
                 cex.axis=1.5, las=2)
            mtext(met1lab, side=2, line=5, cex=1.5)
            if(leg.panel == 'left' & leg.row == 1){
              legend(leg.loc,
                     legend=c('No coevolution, no cospeciation',
                       'No coevolution, cospeciation',
                       'Coevolution, no cospeciation',
                       'Coevolution and cospeciation'),
                     col=topo.col, lty='solid', lwd=1.5, bty='n', cex=1)
            }

          } else if(rule=='quan' & leg.panel == 'right' & leg.row == 1){
            legend(leg.loc,
                   legend=c('No coevolution, no cospeciation',
                     'No coevolution, cospeciation',
                     'Coevolution, no cospeciation',
                     'Coevolution and cospeciation'),
                   col=topo.col, lty='solid', lwd=1.5, bty='n', cex=1)
          }
        }
      })
    }
    plot.met2 <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        if(topos == 'same.same'){
          plot(smooth.spline(met2[,'mean'] ~ intint, df=5), xlab='', ylab='',
               col=topo.col[topo],
               las=1,
               ylim=range.met2,
               xlim=range.intint,
               yaxt='n',
               xaxt='n',
               cex.axis=1.2,
               type='l',
               lwd=2,...)
        }
        not.na <- !is.na(met2[,'ci.ub'])
        out1 <- smooth.spline(intint[not.na], met2[,'ci.lb'][not.na], df=5)
        out2 <- smooth.spline(intint[not.na], met2[,'ci.ub'][not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y

        polygon(x=c(intint[not.na], rev(intint[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)

        points(smooth.spline(met2[,'mean'] ~ intint, df=5),
               col=topo.col[topo],
               type='l', lwd=2, ...)
        if(topos == 'same.same'){
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(rule == 'qual'){
            axis(2, at=pretty(range.met2, n=5, min.n=3),
                 cex.axis=1.5, las=2)
            mtext(met2lab, side=2, line=5, cex=1.5)
            abline(h=0, lty='dashed', col='grey')
            if(leg.row == 2){
              if(leg.panel == 'left'){
                legend(leg.loc,
                       legend=c('Independent evolution',
                         'No coevolution, cospeciation',
                         'Coevolution, no cospeciation',
                         'Coevolution and cospeciation'),
                       col=topo.col, lty='solid', lwd=1.5, bty='n', cex=1)
              }
            }
          }
        }
      })
    }
    plot.met3 <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        if(topos == 'same.same'){
          plot(smooth.spline(met3[,'mean'] ~ intint, df=5), xlab='', ylab='',
               col=topo.col[topo],
               las=1,
               ylim=range.met3,
               xlim=range.intint,
               yaxt='n',
               xaxt='n',
               cex.axis=1.2,
               type='l',
               lwd=2,...)
        }
        not.na <- !is.na(met3[,'ci.ub'])
        out1 <- smooth.spline(intint[not.na], met3[,'ci.lb'][not.na], df=5)
        out2 <- smooth.spline(intint[not.na], met3[,'ci.ub'][not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y
        c(y1.min, rev(y1.max))

        polygon(x=c(intint[not.na], rev(intint[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)

        points(smooth.spline(met3[,'mean'] ~ intint, df=5),
               col=topo.col[topo],
               type='l', lwd=2, ...)

        if(topos == 'same.same'){
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(rule == 'qual'){
            axis(2, at=pretty(range.met3, n=5, min.n=3),
                 cex.axis=1.5, las=2)

            mtext(met3lab, side=2, line=5, cex=1.5)
          }
        }
      })
    }


    layout(matrix(1:6, ncol=2))
    par(oma=c(7, 9, 5, 1), mar=c(0.5, 0.5, 0.5, 1), mgp=c(2, 1, 0))

    plot.met1(res.evo, 'qual', 'same.same')
    plot.met1(res.evo, 'qual', 'diff.same')
    plot.met1(res.evo, 'qual', 'same.diff')
    plot.met1(res.evo, 'qual', 'diff.diff')
    mtext('Unweighted', side=3, line=2, cex=1.5)

    plot.met2(res.evo, 'qual', 'same.same')
    plot.met2(res.evo, 'qual', 'diff.same')
    plot.met2(res.evo, 'qual', 'same.diff')
    plot.met2(res.evo, 'qual', 'diff.diff')

    plot.met3(res.evo, 'qual', 'same.same')
    plot.met3(res.evo, 'qual', 'diff.same')
    plot.met3(res.evo, 'qual', 'same.diff')
    plot.met3(res.evo, 'qual', 'diff.diff')

    axis(1, at=pretty(range.intint, n=5, min.n=3),
         cex.axis=1.5)
    mtext(xlabel, side=1, line=6, cex=1.5)
    mtext(left.lab, side=1, line=4, cex=1.5, adj=0)
    mtext(right.lab, side=1, line=4, cex=1.5, adj=1)

    plot.met1(res.evo, 'quan', 'same.same')
    plot.met1(res.evo, 'quan', 'diff.same')
    plot.met1(res.evo, 'quan', 'same.diff')
    plot.met1(res.evo, 'quan', 'diff.diff')
    mtext('Weighted', side=3, line=2, cex=1.5)

    plot.met2(res.evo, 'quan', 'same.same')
    plot.met2(res.evo, 'quan', 'diff.same')
    plot.met2(res.evo, 'quan', 'same.diff')
    plot.met2(res.evo, 'quan', 'diff.diff')

    plot.met3(res.evo, 'quan', 'same.same')
    plot.met3(res.evo, 'quan', 'diff.same')
    plot.met3(res.evo, 'quan', 'same.diff')
    plot.met3(res.evo, 'quan', 'diff.diff')

    axis(1, at=pretty(range.intint, n=5, min.n=3),
         cex.axis=1.5)
    mtext(xlabel, side=1, line=6, cex=1.5)
    mtext(left.lab, side=1, line=4, cex=1.5, adj=0)
    mtext(right.lab, side=1, line=4, cex=1.5, adj=1)
  }
  pdf.f(f, file= file.path(path.dir, sprintf('%s.pdf',
             paste(strsplit(xmetric, '\\.')[[1]][1],
                   metric1, metric2, metric3,
                   sep=''))), width=8, height=9)

}
