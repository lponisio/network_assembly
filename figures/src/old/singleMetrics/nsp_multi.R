library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
nsp.plot <- function(simres, metric1, metric2, pp, ylabs){
  f <- function(){
    res.sim <- aggregate(list(NODF= simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           nsp = simres[[1]][, pp]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.nodf <- aggregate(list(CI=simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           nsp = simres[[1]][, pp]), 
                         function(x) sd(x, na.rm=TRUE))
    
    res.sim$nodf.ci.lb <-   res.sim$NODF + sd.nodf$CI
    res.sim$nodf.ci.ub <-   res.sim$NODF - sd.nodf$CI

    mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                     by= list(topo= simres[[2]][,"topo"],
                       link.rule= simres[[2]][,"link.rule"],
                       mat= simres[[2]][,"mats"],
                       nsp = simres[[1]][, pp]), 
                     function(x) mean(x, na.rm=TRUE))

    sd.mod <- aggregate(list(CI=simres[[1]][, metric2]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"],
                          nsp = simres[[1]][, pp]), 
                        function(x) sd(x, na.rm=TRUE)
                        )
    res.sim$mod <- mod$mod
    res.sim$mod.ci.lb <-   res.sim$mod + sd.mod$CI
    res.sim$mod.ci.ub <-   res.sim$mod - sd.mod$CI


    range.NODF <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub),
                        na.rm=TRUE)
    range.mod <- range(c(res.sim$mod.ci.lb, res.sim$mod.ci.ub), na.rm=TRUE)
    range.nsp <- range(res.sim$nsp, na.rm =TRUE)
    
    res.evo <- res.sim[res.sim$mat == 'evo',]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')

    layout(matrix(1:6, ncol=3, nrow=2, byrow= TRUE))
    
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(6,8,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2, yaxt="n")
      axis(2, at=pretty(range.NODF, n=3), cex.axis=1.5, las=2)
      ## arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
      mtext('Matching', side=3, line=2, cex=1.5)
      mtext('Nestedness', side=2, line=5, cex=1.5)
    })
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2, yaxt="n")
      ## arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
      mtext('Barrier', side=3, line=2, cex=1.5)
    })
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(nsp, NODF, ylim=range.NODF, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp,
           las=2, xaxt="n", yaxt="n")
      ## arrows(y0=nodf.ci.lb, y1=nodf.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
      mtext('Neutral', side=3, line=2, cex=1.5)
    })
    legend("topright",
           legend=c('Independent evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bty="n", cex=1.5)
    
    with(res.evo[res.evo$link.rule == 'matching',],{
      plot(nsp, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2, yaxt="n")
      axis(1, at= pretty(range.nsp, n=3),  cex.axis=1.5)
      ## arrows(y0=mod.ci.lb, y1=mod.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
      axis(2, at=pretty(range.mod, n=4), cex.axis=1.5, las=2)
      mtext('Modularity', side=2, line=5, cex=1.5)
    })
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(nsp, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp, xaxt='n',
           las=2, yaxt="n")
      axis(1, at= pretty(range.nsp, n=3),  cex.axis=1.5)
      ## arrows(y0=mod.ci.lb, y1=mod.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
    })
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(nsp, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.nsp,
           las=2, xaxt="n", yaxt="n")
      axis(1, at= pretty(range.nsp, n=3),  cex.axis=1.5)
      ## arrows(y0=mod.ci.lb, y1=mod.ci.ub, x0=nsp, angle=90,
      ##        length=0, code=0,  col=topo.col[topo])
    })
    
    mtext('Number of species', side=1, line=4, outer=TRUE, cex=1.5)
  }
  path <- "~/Dropbox/network_assembly/figures/coal/singleMetrics"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(pp, metric1, metric2,
             sep=""))), width=9, height=6)
  
}

