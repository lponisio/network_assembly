intInt.plot2 <- function(simres,
                         metric1,
                         metric2,
                         met1lab='Nestedness',
                         met2lab='Modularity',
                         xmetric,
                         xlabel=NA,
                         path.dir,
                         subset=FALSE,
                         column=NA,
                         case=NA,
                         leg.panel='left',
                         leg.loc='bottomleft',
                         drop.same.same=FALSE,
                         left.lab='Low',
                         right.lab='High'){
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
                            met2= simres[[1]][, metric2]),
                       by= list(topo= simres[[2]][,'topo'],
                         link.rule= simres[[2]][,'link.rule'],
                         mat= simres[[2]][,'mats'],
                         intint = simres[[1]][, xmetric]), 
                       calc.stats)
  res.sim <- res.sim[!is.na(res.sim$met1[,'ci.lb']),]
  res.sim <- res.sim[!is.na(res.sim$met2[,'ci.lb']),]
  
  ## metric 1
  range.met1 <- range(c(res.sim$met1[,'ci.lb'],
                        res.sim$met1[,'ci.ub']),
                      na.rm=TRUE)
  ## metric 2
  range.met2 <- range(c(res.sim$met2[,'ci.lb'],
                        res.sim$met2[,'ci.ub']),
                      na.rm=TRUE)

  res.evo <- res.sim[res.sim$mat == 'evo',]
  res.evo <- res.evo[!apply(res.evo, 1, function(x) any(is.na(x))),]
  
  ## if(xmetric == 'range.size'){
  ##   range.intint <- rev(range(res.evo$intint, na.rm =TRUE))
  ## } else {
    range.intint <- range(res.evo$intint, na.rm =TRUE)
  
  topo.col <- brewer.pal(11, 'Spectral')[c(1,3,10,11)]
  fill.col <- add.alpha(topo.col, alpha=0.2)
  names(topo.col) <- names(fill.col) <- 
    c('diff.diff','same.diff','diff.same','same.same')
  f <- function(){
    layout(matrix(1:4, ncol=2))
    par(oma=c(7, 9, 5, 1), mar=c(0.5, 0.5, 0.5, 1), mgp=c(2, 1, 0))
    plot.met1 <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        if(topos == 'diff.same'){
          plot(smooth.spline(met1[,'mean'] ~ intint, df=5),
               xlab='', ylab='',
               col=topo.col[topo],
               las=1,
               ylim=range.met1,
               xlim=range.intint,
               yaxt='n',
               xaxt='n',
               type='l',
               lwd=2,...)
        }
        not.na <- !is.na(met1[,'ci.ub'])
        out1 <- smooth.spline(intint[not.na], met1[,'ci.lb'][not.na], df=5)
        out2 <- smooth.spline(intint[not.na], met1[,'ci.ub'][not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y

        polygon(x=c(intint[not.na], rev(intint[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)
        points(smooth.spline(met1[,'mean'] ~ intint, df=5),
               col=topo.col[topo], 
               type='l', lwd=2, ...)
        if(topos == 'diff.same'){
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(rule == 'qual'){
            axis(2, at=pretty(range.met1, n=5, min.n=3),
                 cex.axis=1.5, las=2)
            mtext(met1lab, side=2, line=5, cex=1.5)
          }
        }
      })
    }

    plot.met2 <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        if(topos == 'diff.same'){
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
        
        if(topos == 'diff.same'){
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(rule == 'qual'){
            axis(2, at=pretty(range.met2, n=5, min.n=3),
                 cex.axis=1.5, las=2)
            if(leg.panel == 'left'){
              legend(leg.loc,
                     legend=c('No coevolution, no cospeciation',
                       'No coevolution, cospeciation',
                       'Coevolution, no cospeciation',
                       'Coevolution and cospeciation'),
                     col=topo.col, lty='solid', lwd=1.5, bty='n', cex=0.6)
            }
            mtext(met2lab, side=2, line=5, cex=1.5)
          } else {
            if(leg.panel == 'right'){
              legend(leg.loc,
                     legend=c('No coevolution, no cospeciation',
                       'No coevolution, cospeciation',
                       'Coevolution, no cospeciation',
                       'Coevolution and cospeciation'),
                     col=topo.col, lty='solid', lwd=1.5, bty='n', cex=0.6)
            }
          }
        }
      })
    }
    
    plot.met1(res.evo, 'qual', 'diff.same')
    plot.met1(res.evo, 'qual', 'same.diff')
    plot.met1(res.evo, 'qual', 'diff.diff')
    if(!drop.same.same)  plot.met1(res.evo, 'qual', 'same.same')
    mtext('Unweighted', side=3, line=2, cex=1.5)

    plot.met2(res.evo, 'qual', 'diff.same')
    plot.met2(res.evo, 'qual', 'same.diff')
    plot.met2(res.evo, 'qual', 'diff.diff')
    if(!drop.same.same) plot.met2(res.evo, 'qual', 'same.same')
    
    if(xmetric=='range.size'){
      mtext(left.lab, side=1, line=3, cex=1.2, adj=0)
      mtext(right.lab, side=1, line=3, cex=1.2, adj=1)
      axis(side=1, pretty(range.intint, 4), cex.axis=1.5)
    } else{
      axis(side=1, pretty(range.intint, 4), cex.axis=1.5)
    }
    mtext(xlabel, side=1, line=5, cex=1.5)

    plot.met1(res.evo, 'quan', 'diff.same')
    plot.met1(res.evo, 'quan', 'same.diff')
    plot.met1(res.evo, 'quan', 'diff.diff')
    if(!drop.same.same) plot.met1(res.evo, 'quan', 'same.same')
    mtext('Weighted', side=3, line=2, cex=1.5)

    plot.met2(res.evo, 'quan', 'diff.same')
    plot.met2(res.evo, 'quan', 'same.diff')
    plot.met2(res.evo, 'quan', 'diff.diff')
    if(!drop.same.same) plot.met2(res.evo, 'quan', 'same.same')

    mtext(xlabel, side=1, line=5, cex=1.5)## , adj=-0.7
    
    if(xmetric=='range.size'){
      mtext(left.lab, side=1, line=3, cex=1.2, adj=0)
      mtext(right.lab, side=1, line=3, cex=1.2, adj=1)
      axis(side=1, pretty(range.intint, 4), cex.axis=1.5)
    } else{
      axis(side=1, pretty(range.intint, 4), cex.axis=1.5)
    }
  }
  pdf.f(f, file= file.path(path.dir, sprintf('%s.pdf',
             paste(strsplit(xmetric, '\\.')[[1]][1],
                   metric1,
                   metric2, sep=''))),
        width=7, height=6)
  
}
