library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, fig.height, path, metric1,
                             metric2){

  calc.overlap <- function(ranges, res.sim){
    mins.nodf1 <- matrix(rep(ranges[,1], each =
                             nrow(res.sim)), nrow=nrow(res.sim))
    mins.nodf2 <- t(mins.nodf1)
    mins.tot <- pmax(mins.nodf1, mins.nodf2)
    maxs.nodf1 <- matrix(rep(ranges[,2], each =
                             nrow(res.sim)), nrow=nrow(res.sim))
    maxs.nodf2 <- t(maxs.nodf1)
    maxs.tot <- pmin(maxs.nodf1, maxs.nodf2)
    overlap <- maxs.tot - mins.tot
    colnames(overlap) <- res.sim$topo
    rownames(overlap) <- res.sim$topo
    return(overlap)
  }
  f <- function(){
    
    res.sim <- aggregate(list(NODF= simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.nodf <- aggregate(list(CI=simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) range(x, na.rm=TRUE))
    res.sim$nodf.ci.lb <-  sd.nodf$CI[,1]
    res.sim$nodf.ci.ub <-  sd.nodf$CI[,2]
    means.mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) mean(x, na.rm=TRUE))
    sd.mod <- aggregate(list(CI= simres[[1]][, metric2]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"]), 
                        function(x) range(x, na.rm=TRUE))
    res.sim$mod <- means.mod$mod
    res.sim$mod.ci.lb <-  sd.mod$CI[,1]
    res.sim$mod.ci.ub <-  sd.mod$CI[,2]
    
    ## range.NODF <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub))
    ## range.mod <- range(c(res.sim$mod.ci.lb, res.sim$mod.ci.ub)) 
    ## evo mats
    res.evo <- res.sim[res.sim$mat == 'evo',]
    res.evo.b <- res.evo[res.evo$link.rule == 'barrior',]
    res.evo.m <- res.evo[res.evo$link.rule == 'matching',]
    res.evo.n <- res.evo[res.evo$link.rule == 'neutral',]
    overlap.evo.b.nodf <- calc.overlap(cbind(res.evo.b$nodf.ci.lb,
                                             res.evo.b$nodf.ci.ub),
                                       res.sim= res.evo.b)
    overlap.evo.b.mod <- calc.overlap(cbind(res.evo.b$mod.ci.lb,
                                            res.evo.b$mod.ci.ub),
                                      res.sim= res.evo.b)
    overlap.evo.m.mod <- calc.overlap(cbind(res.evo.m$mod.ci.lb,
                                            res.evo.m$mod.ci.ub),
                                      res.sim= res.evo.m)
    overlap.evo.m.nodf <- calc.overlap(cbind(res.evo.m$nodf.ci.lb,
                                             res.evo.m$nodf.ci.ub),
                                       res.sim= res.evo.m)
    overlap.evo.n.nodf <- calc.overlap(cbind(res.evo.n$nodf.ci.lb,
                                             res.evo.n$nodf.ci.ub),
                                       res.sim= res.evo.n)
    overlap.evo.n.mod <- calc.overlap(cbind(res.evo.n$mod.ci.lb,
                                            res.evo.n$mod.ci.ub),
                                      res.sim= res.evo.n)
    
    res.eco <- res.sim[res.sim$mat == 'eco',]
    res.eco.b <- res.eco[res.eco$link.rule == 'barrior',]
    res.eco.m <- res.eco[res.eco$link.rule == 'matching',]
    res.eco.n <- res.eco[res.eco$link.rule == 'neutral',]
    overlap.eco.b.nodf <- calc.overlap(cbind(res.eco.b$nodf.ci.lb,
                                             res.eco.b$nodf.ci.ub),
                                       res.sim= res.eco.b)
    overlap.eco.b.mod <- calc.overlap(cbind(res.eco.b$mod.ci.lb,
                                            res.eco.b$mod.ci.ub),
                                      res.sim= res.eco.b)
    overlap.eco.m.mod <- calc.overlap(cbind(res.eco.m$mod.ci.lb,
                                            res.eco.m$mod.ci.ub),
                                      res.sim= res.eco.m)
    overlap.eco.m.nodf <- calc.overlap(cbind(res.eco.m$nodf.ci.lb,
                                             res.eco.m$nodf.ci.ub),
                                       res.sim= res.eco.m)
    overlap.eco.n.nodf <- calc.overlap(cbind(res.eco.n$nodf.ci.lb,
                                             res.eco.n$nodf.ci.ub),
                                       res.sim= res.eco.n)
    overlap.eco.n.mod <- calc.overlap(cbind(res.eco.n$mod.ci.lb,
                                            res.eco.n$mod.ci.ub),
                                      res.sim= res.eco.n)
    
    res.eco2 <- res.sim[res.sim$mat == 'eco2',]
    res.eco2.b <- res.eco2[res.eco2$link.rule == 'barrior',]
    res.eco2.m <- res.eco2[res.eco2$link.rule == 'matching',]
    res.eco2.n <- res.eco2[res.eco2$link.rule == 'neutral',]
    overlap.eco2.b.nodf <- calc.overlap(cbind(res.eco2.b$nodf.ci.lb,
                                              res.eco2.b$nodf.ci.ub),
                                        res.sim= res.eco2.b)
    overlap.eco2.b.mod <- calc.overlap(cbind(res.eco2.b$mod.ci.lb,
                                             res.eco2.b$mod.ci.ub),
                                       res.sim= res.eco2.b)
    overlap.eco2.m.mod <- calc.overlap(cbind(res.eco2.m$mod.ci.lb,
                                             res.eco2.m$mod.ci.ub),
                                       res.sim= res.eco2.m)
    overlap.eco2.m.nodf <- calc.overlap(cbind(res.eco2.m$nodf.ci.lb,
                                              res.eco2.m$nodf.ci.ub),
                                        res.sim= res.eco2.m)
    overlap.eco2.n.nodf <- calc.overlap(cbind(res.eco2.n$nodf.ci.lb,
                                              res.eco2.n$nodf.ci.ub),
                                        res.sim= res.eco2.n)
    overlap.eco2.n.mod <- calc.overlap(cbind(res.eco2.n$mod.ci.lb,
                                             res.eco2.n$mod.ci.ub),
                                       res.sim= res.eco2.n)
    
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <- c('diff.diff','same.diff','diff.same','same.same')

    layout(matrix(1:6, ncol=2, nrow=3, byrow= TRUE))
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n',
           las=2)
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Matching', side=2, line=5)
      mtext('Evolutionary', side=3, line=2)
    })
    with(res.eco[res.eco$link.rule == 'matching',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab="", yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n')
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Ecological', side=3, line=2)
    }) 
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF,
           xaxt='n', las=2)
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Barrier', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'barrior',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n')
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
    })
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, las=1)
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Neutral', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'neutral',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF)
      arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
             angle=90, length=0.05, code=3)
      arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      legend("topright",
             legend=c('Neutral evolution','No coevolution, co-speciation',
               'Coevolution, no co-speciation','Coevolution and co-speciation'),
             col=topo.col, pch=16, bty="n")
    })

    mtext('Nestedness', side=1, line=3, outer=TRUE)
    mtext('Modularity', side=2, line=3, outer=TRUE)
  }

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1, metric2,
             sep=""))), width=7, height=7)
  
}

