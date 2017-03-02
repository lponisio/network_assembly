library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
corInt.plot <- function(simres, fig.height, path){
  f <- function(){
    res.sim <- aggregate(list(cor.pol= simres[[1]][, "cor.pol"]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.pol <- aggregate(list(cor.pol= simres[[1]][, "cor.pol"]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"]), 
                        function(x) sd(x, na.rm=TRUE))
    res.sim$sd.pol <- sd.pol$cor.pol
    means.plant <- aggregate(list(cor.plant= simres[[1]][, "cor.plant"]), 
                             by= list(topo= simres[[2]][,"topo"],
                               link.rule= simres[[2]][,"link.rule"],
                               mat= simres[[2]][,"mats"]), 
                             function(x) mean(x, na.rm=TRUE))
    sd.plant <- aggregate(list(cor.plant= simres[[1]][, "cor.plant"]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) sd(x, na.rm=TRUE))
    res.sim$cor.plant <- means.plant$cor.plant
    res.sim$sd.plant <- sd.plant$cor.plant
    res.sim <- res.sim[is.finite(res.sim$cor.pol),]
    res.evo <- res.sim[res.sim$mat == 'evo',]
    res.eco <- res.sim[res.sim$mat == 'eco',]
    res.eco2 <- res.sim[res.sim$mat == 'eco2',]

    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <- c('diff.diff','same.diff','diff.same','same.same')
    layout(matrix(1:9, ncol=3, nrow=3, byrow= TRUE))


    ##**********************************************************
    ## matching
    ##**********************************************************
    
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), xaxt='n',
           las=2)
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90, length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Matching', side=2, line=5)
      mtext('Evolutionary', side=3, line=2.5)
    })
    with(res.eco[res.eco$link.rule == 'matching',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab="",
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), xaxt='n')
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90, length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Evo-Evological', side=3, line=2.5)
      mtext('Random Abundances', side=3, line=1)
    })

    with(res.eco2[res.eco2$link.rule == 'matching',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab="",
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), xaxt='n')
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Evo-Evological', side=3, line=2.5)
      mtext('Phylo Abundances', side=3, line=1)
    })

    
    ##**********************************************************
    ## barrier
    ##**********************************************************
    
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8),
           xaxt='n', las=2)
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Barrier', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'barrior',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), xaxt='n')
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
    })

    with(res.eco2[res.eco2$link.rule == 'barrior',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), xaxt='n')
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
    })


    ##**********************************************************
    ## neutral
    ##**********************************************************  
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8), las=1)
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Neutral', side=2, line=5)
      legend("topleft",
             legend=c('Neutral evolution','No coevolution, co-speciation',
               'Coevolution, no co-speciation',
      'Coevolution and co-speciation'), col=topo.col, pch=16, bty="n")
    })
    with(res.eco[res.eco$link.rule == 'neutral',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8))
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
    })

    with(res.eco2[res.eco2$link.rule == 'neutral',],{
      plot(cor.plant, cor.pol, ylim=c(-0.1,0.8), xlab='', ylab='',
           yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xlim=c(-0.1,0.8))
      arrows(x0=cor.plant, y0=cor.pol-sd.pol, y1=cor.pol+sd.pol,
             angle=90, length=0.05, code=3)
      arrows(x0=cor.plant-sd.plant, x1=cor.plant+sd.plant, y0=cor.pol,
             angle=90,
             length=0.05, code=3)
      points(cor.plant, cor.pol,
             col=topo.col[topo], pch=16, cex=2)
    })

    mtext('Plant/host phylogenetic interaction signal', side=1,
          line=3, outer=TRUE)
    mtext('Animal/parasite phylogenetic interaction signal', side=2,
          line=3, outer=TRUE)
    
  }

  pdf.f(f, file= file.path(path, "corInt.pdf"), width=7, height=7)
  
}

