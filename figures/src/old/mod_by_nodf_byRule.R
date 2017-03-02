library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, fig.height, path, metric1, metric2){
  f <- function(){

    res.sim <- aggregate(list(NODF= simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.nodf <- aggregate(list(NODF= simres[[1]][, metric1]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) sd(x, na.rm=TRUE))
    res.sim$sd.NODF <- sd.nodf$NODF
    means.mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) mean(x, na.rm=TRUE))
    sd.mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"]), 
                        function(x) sd(x, na.rm=TRUE))
    res.sim$mod <- means.mod$mod
    res.sim$sd.mod <- sd.mod$mod
    res.evo <- res.sim[res.sim$mat == 'evo',]
    res.eco <- res.sim[res.sim$mat == 'eco',]

    range.NODF <- range(simres[[1]][, metric1], na.rm=TRUE)
    range.mod <- range(simres[[1]][, metric2], na.rm=TRUE)

    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <- c('diff.diff','same.diff','diff.same','same.same')
    layout(matrix(1:6, ncol=2, nrow=3, byrow= TRUE))
    
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n',
           las=2)
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Matching', side=2, line=5)
      mtext('Evolutionary', side=3, line=2)
    })
    with(res.eco[res.eco$link.rule == 'matching',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab="", yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n')
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Ecological', side=3, line=2)
    })
    
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF,
           xaxt='n', las=2)
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Barrier', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'barrior',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, xaxt='n')
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
    })
    
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(NODF, mod, ylim=range.mod, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF, las=1)
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3)
      points(NODF, mod,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Neutral', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'neutral',],{
      plot(NODF,mod, ylim=range.mod, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=range.NODF)
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
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
  path <- "~/Dropbox/network_assembly/figures/coal/NODFbyMod"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1, metric2,
             sep=""))), width=7, height=7)
  
}

