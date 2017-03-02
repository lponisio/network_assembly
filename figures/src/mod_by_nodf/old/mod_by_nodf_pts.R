library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, fig.height, path, metric1, metric2){
  f <- function(){
    sim.res <- as.data.frame(simres[[1]])
    sim.res <- cbind(sim.res, simres[[2]])

    colnames(sim.res)[colnames(sim.res) == metric1] <- "metric1"
    colnames(sim.res)[colnames(sim.res) == metric2] <- "metric2"
    
    range.NODF <- quantile(sim.res[,"metric1"], probs= c(0.01,
                                                  0.99), na.rm=TRUE)
    range.mod <- quantile(sim.res[,"metric2"], probs= c(0.01,
                                                 0.99),
                                                  na.rm=TRUE)
    ## range.mod <- range(sim.res[,"metric2"], na.rm=TRUE)
    ## range.NODF <- range(sim.res[,"metric1"], na.rm=TRUE)
    
    res.evo <- sim.res[sim.res$mats == 'evo',]
    res.eco <- sim.res[sim.res$mats == 'eco',]
    res.eco2 <- sim.res[sim.res$mats == 'eco2',]
    
    ## topo.col <- brewer.pal(4, 'Spectral')
    topo.col <-hsv(c(0, 0.2, 0.4,0.6), alpha=0.1)
    names(topo.col) <- c('diff.diff','same.diff','diff.same','same.same')

    ##**********************************************************
    ## matching
    ##**********************************************************
    
    layout(matrix(1:9, ncol=3, nrow=3, byrow= TRUE))
    with(res.evo[res.evo$link.rule == 'matching',],{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
      plot(metric1, metric2, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt='n',
           las=2, ylim=range.mod, xlim=range.NODF)
      mtext('Matching', side=2, line=5)
      mtext('Evolutionary', side=3, line=2.2)
    })
    with(res.eco[res.eco$link.rule == 'matching',],{
      plot(metric1, metric2, xlab='', ylab="", yaxt="n",
           col=topo.col[topo], pch=16, cex=2, xaxt='n',
           xlim=range.NODF, ylim=range.mod)
      mtext('Evo-Evological', side=3, line=2.5)
      mtext('Random Abundances', side=3, line=1)
    })
    with(res.eco2[res.eco2$link.rule == 'matching',],{
      plot(metric1, metric2 ,xlab='', ylab="", yaxt="n",
           col=topo.col[topo], pch=16, cex=2,
           xaxt='n', ylim=range.mod, xlim=range.NODF)
      mtext('Evo-Evological', side=3, line=2.5)
      mtext('Phylo Abundances', side=3, line=1)
    })

    ##**********************************************************
    ## barrier
    ##**********************************************************
    
    with(res.evo[res.evo$link.rule == 'barrior',],{
      plot(metric1, metric2, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2,
           xaxt='n', las=2, ylim=range.mod,  xlim=range.NODF)
      mtext('Barrier', side=2, line=5)
    })
    with(res.eco[res.eco$link.rule == 'barrior',],{
      plot(metric1, metric2, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2, 
           xaxt='n', ylim=range.mod, xlim=range.NODF)
    })
    with(res.eco2[res.eco2$link.rule == 'barrior',],{
      plot(metric1, metric2, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2,
           xaxt='n', ylim=range.mod,  xlim=range.NODF)
    })

    ##**********************************************************
    ## neutral
    ##**********************************************************
    
    with(res.evo[res.evo$link.rule == 'neutral',],{
      plot(metric1, metric2, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, las=1,
           ylim=range.mod, xlim=range.NODF)
      mtext('Neutral', side=2, line=5)
      legend("topleft",
             legend=c('Neutral evolution','No coevolution, co-speciation',
               'Coevolution, no co-speciation',
               'Coevolution and co-speciation'),
             col=topo.col, pch=16, bty="n")
    })
    with(res.eco[res.eco$link.rule == 'neutral',],{
      plot(metric1, metric2, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2,
           ylim=range.mod, xlim=range.NODF)
    })
    with(res.eco2[res.eco2$link.rule == 'neutral',],{
      plot(metric1, metric2, xlab='', ylab='', yaxt="n",
           col=topo.col[topo], pch=16, cex=2,
           xlim=range.NODF,  ylim=range.mod)
    })
    mtext('Nestedness', side=1, line=3, outer=TRUE)
    mtext('Modularity', side=2, line=3, outer=TRUE)
  }

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1, metric2,"pts",
             sep=""))), width=7, height=7)
  
}

