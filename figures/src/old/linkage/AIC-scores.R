library(RColorBrewer)

aic.plot <- function(simres, metric){

  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  f <- function(){
    simres <- simres[, c("topo", "mats", "link.rule", metric)]
    lr <- numeric(nrow(simres))
    for(i in 1:nrow(simres)){
      lr[i] <-  ifelse(simres[i,"link.rule"] == "barrior", 1,2)
    }

    simres$lr <- lr
    colnames(simres)[colnames(simres) == metric] <- "metric"
    
    res.evo <- simres[simres$mats == 'evo',]
    res.eco <- simres[simres$mats == 'eco',]
    res.eco2 <- simres[simres$mats == 'eco2',]
    
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    
    layout(matrix(1:3, ncol=3, nrow=1, byrow= TRUE))
    with(res.evo,{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,2),mgp=c(2,1,0))
      plot(x=lr, y=metric, ylim=c(0,1), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt="n",
           las=2, xlim=c(0,3))
      mtext('Evolutionary', side=3, line=2.5)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
                                        #mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
                                        #abline(v=2.5, lty="dashed")
    })
    with(res.eco,{
      plot(x=lr, y=metric, ylim=c(0,1), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt='n', yaxt="n",
           xlim=c(0,3))
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Random abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
                                        #mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
                                        #abline(v=2.5, lty="dashed")
    })
    with(res.eco2,{
      plot(x=lr, y=metric, ylim=c(0,1), xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt='n', yaxt="n",
           xlim=c(0,3))              
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Phylo abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
                                        #mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
                                        # abline(v=2.5, lty="dashed")
      ## legend("topright",
      ##        legend=c('Neutral evolution','No coevolution, co-speciation',
      ##          'Coevolution, no co-speciation','Coevolution and co-speciation'),
      ##        col=topo.col, pch=16, bty="n")
    })
  }
  path <- "~/Dropbox/network_assembly/figures/coal/SingleMetrics"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("aic", metric,
             sep="_"))), width=9, height=3)

}

