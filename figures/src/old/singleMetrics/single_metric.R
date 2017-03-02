library(RColorBrewer)

single.met.plot <- function(simres, metric){

  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  f <- function(){
    res.sim <- aggregate(list(metric= simres[[1]][, metric]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    
    CI.metric <- aggregate(list(CI=simres[[1]][, metric]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) quantile(x, probs= c(0.025,
                                                     0.975), na.rm=TRUE))
    res.sim$ci.lb <-  CI.metric$CI[,1]
    res.sim$ci.ub <-  CI.metric$CI[,2]
    res.sim$lr <-rep(rep(1:3, each=4),3)
    
    range.metric <- range(c(res.sim$ci.lb, res.sim$ci.ub), na.rm=TRUE)
    
    res.evo <- res.sim[res.sim$mat == 'evo',]
    res.eco <- res.sim[res.sim$mat == 'eco',]
    res.eco2 <- res.sim[res.sim$mat == 'eco2',]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    layout(matrix(1:3, ncol=3, nrow=1, byrow= TRUE))
    with(res.evo,{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,2),mgp=c(2,1,0))
      plot(x=lr, y=metric, ylim=range.metric, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt="n",
           las=2)
      arrows(y0=ci.lb, y1=ci.ub, x0=lr, angle=90,
             length=0.05, code=3)
      points(lr, metric,
             col=topo.col[topo], pch=16, cex=2)
                                       
      mtext('Evolutionary', side=3, line=2.5)
      mtext('Matching', side=1, line=2, adj=0)
      mtext('Barrier', side=1, line=2)
      mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
      abline(v=2.5, lty="dashed")
    })
    with(res.eco,{
      plot(x=lr, y=metric, ylim=range.metric, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt='n', yaxt="n")
      arrows(y0=ci.lb, y1=ci.ub, x0=lr, angle=90,
             length=0.05, code=3)
      points(lr, metric,
             col=topo.col[topo], pch=16, cex=2)
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Random abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0)
      mtext('Barrier', side=1, line=2)
      mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
      abline(v=2.5, lty="dashed")
    })
    with(res.eco2,{
      plot(x=lr, y=metric, ylim=range.metric, xlab='', ylab='',
           col=topo.col[topo], pch=16, cex=2, xaxt='n', yaxt="n")
      arrows(y0=ci.lb, y1=ci.ub, x0=lr, angle=90,
             length=0.05, code=3)
      points(lr, metric,
             col=topo.col[topo], pch=16, cex=2)               
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Phylo abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0)
      mtext('Barrier', side=1, line=2)
      mtext('Neutral', side=1, line=2, adj=1)
      abline(v=1.5, lty="dashed")
      abline(v=2.5, lty="dashed")
      ## legend("topright",
      ##        legend=c('Neutral evolution','No coevolution, co-speciation',
      ##          'Coevolution, no co-speciation','Coevolution and co-speciation'),
      ##        col=topo.col, pch=16, bty="n")
    })
  }
  path <- "~/Dropbox/network_assembly/figures/coal/SingleMetrics"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric,
             sep=""))), width=9, height=3)

}

