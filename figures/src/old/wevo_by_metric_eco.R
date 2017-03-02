library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
w.evo.plot <- function(simres, fig.height, metric){
  f <- function(){
    res.sim <- aggregate(list(NODF= simres[[1]][, metric]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           w.evo = simres[[1]][, "w.evo"]), 
                         function(x) mean(x, na.rm=TRUE))
    
    sd.nodf <- aggregate(list(CI=simres[[1]][, metric]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           w.evo = simres[[1]][, "w.evo"]), 
                         function(x) quantile(x, probs= c(0.025,
                                                   0.975), na.rm=TRUE))
    res.sim$nodf.ci.lb <-  sd.nodf$CI[,1]
    res.sim$nodf.ci.ub <-  sd.nodf$CI[,2]
    
    range.NODF <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub), na.rm=TRUE)
    range.wevo <- range(res.sim$w.evo, na.rm =TRUE)
    
    res.evo <- res.sim[res.sim$mat == 'eco',]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')

    layout(matrix(1:3, ncol=1, nrow=3, byrow= TRUE))
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(w.evo, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.wevo, xaxt='n',
           las=2, log='x')
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=w.evo, angle=90,
             length=0.05, code=3)
      points(w.evo, NODF,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Matching', side=2, line=5)
      mtext('Ecological', side=3, line=2)
    })
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(w.evo, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.wevo, xaxt='n',
           las=2,  log='x')
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=w.evo, angle=90,
             length=0.05, code=3)
      points(w.evo, NODF,
             col=topo.col[topo], pch=16, cex=2)
       mtext('Barrier', side=2, line=5)
    })
    with(res.evo[res.eco$link.rule == 'neutral',],{
      plot(w.evo, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.wevo,
           las=2,  log='x')
      arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=w.evo, angle=90,
             length=0.05, code=3)
      points(w.evo, NODF,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Neutral', side=2, line=5)
      legend("topright",
             legend=c('Neutral evolution','No coevolution, co-speciation',
               'Coevolution, no co-speciation','Coevolution and co-speciation'),
             col=topo.col, pch=16, bty="n")
    })
    
    mtext('Evolutionary weight', side=1, line=3, outer=TRUE)
    mtext(metric, side=2, line=3, outer=TRUE)
  }
  path <- "~/Dropbox/network_assembly/figures/coal/wevo"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("wevo", metric,
             sep=""))), width=7, height=7)
  
}

