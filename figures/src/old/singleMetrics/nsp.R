library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
nsp.plot <- function(simres, metric, pp, ylabs){
  f <- function(){
    res.sim <- aggregate(list(NODF= simres[[1]][, metric]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           nsp = simres[[1]][, pp]), 
                         function(x) mean(x, na.rm=TRUE))

    sd.nodf <- aggregate(list(CI=simres[[1]][, metric]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           nsp = simres[[1]][, pp]), 
                         function(x) sd(x, na.rm=TRUE)
                         )
    
    res.sim$nodf.ci.lb <-   res.sim$NODF + sd.nodf$CI
    res.sim$nodf.ci.ub <-   res.sim$NODF - sd.nodf$CI
    
    range.NODF <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub), na.rm=TRUE)
    range.nsp <- range(res.sim$nsp, na.rm =TRUE)
    
    res.evo <- res.sim[res.sim$mat == 'evo',]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')

    layout(matrix(1:3, ncol=1, nrow=3, byrow= TRUE))
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2)
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
              length=0, code=0,  col=topo.col[topo])
      mtext('Matching', side=2, line=5)
     ## mtext('No species abundance', side=3, line=2)
    })
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2)
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
             length=0, code=0,  col=topo.col[topo])
      mtext('Barrier', side=2, line=5)
    })
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp,
           las=2, xaxt="n")
      axis(1, at=seq(from=15,to=105,length=7))
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
              length=0, code=0,  col=topo.col[topo])
      mtext('Neutral', side=2, line=5)
      legend("topright",
             legend=c('Neutral evolution','No coevolution, cospeciation',
               'Coevolution, no cospeciation','Coevolution and cospeciation'),
             col=topo.col, pch=16, bty="n", cex=1.5)
    })
    
    mtext('Number of species', side=1, line=3, outer=TRUE)
    mtext(ylabs, side=2, line=3, outer=TRUE)
  }
  path <- "~/Dropbox/network_assembly/figures/coal/singleMetrics"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(pp, metric,
             sep=""))), width=5, height=7)
  
}

