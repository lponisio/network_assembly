library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mods.plot <- function(simres, path, metric1,
                      metric2, metric3){
  res.sim <- aggregate(list(metric1= simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) mean(x, na.rm=TRUE))
  sd.metric1 <- aggregate(list(CI=simres[[1]][, metric1]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) sd(x, na.rm=TRUE)
                          )
  res.sim$metric1.ci.lb <-   res.sim$metric1 + sd.metric1$CI
  res.sim$metric1.ci.ub <-   res.sim$metric1 - sd.metric1$CI
  
  means.metric2 <- aggregate(list(metric2= simres[[1]][, metric2]), 
                             by= list(topo= simres[[2]][,"topo"],
                               link.rule= simres[[2]][,"link.rule"],
                               mat= simres[[2]][,"mats"]), 
                             function(x) mean(x, na.rm=TRUE))
  sd.metric2 <- aggregate(list(CI= simres[[1]][, metric2]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) sd(x, na.rm=TRUE)
                          )
  res.sim$metric2 <- means.metric2$metric2
  res.sim$metric2.ci.lb <-   res.sim$metric2 + sd.metric2$CI 
  res.sim$metric2.ci.ub <-   res.sim$metric2 - sd.metric2$CI
  
  means.metric3 <- aggregate(list(metric3= simres[[1]][, metric3]), 
                             by= list(topo= simres[[2]][,"topo"],
                               link.rule= simres[[2]][,"link.rule"],
                               mat= simres[[2]][,"mats"]), 
                             function(x) mean(x, na.rm=TRUE))
  sd.metric3 <- aggregate(list(CI= simres[[1]][, metric3]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) sd(x, na.rm=TRUE)
                          )
  res.sim$metric3 <- means.metric3$metric3
  res.sim$metric3.ci.lb <-   res.sim$metric3 + sd.metric3$CI 
  res.sim$metric3.ci.ub <-   res.sim$metric3 - sd.metric3$CI
  
  range.metrics <- range(c(res.sim$metric1.ci.lb,
                           res.sim$metric2.ci.lb,
                           res.sim$metric3.ci.lb,
                           res.sim$metric1.ci.ub,
                           res.sim$metric2.ci.ub,
                           res.sim$metric3.ci.ub))

  res.evo <- res.sim[res.sim$mat == 'evo',]
  topo.col <- brewer.pal(4, 'Spectral')
  names(topo.col) <-
    c('diff.diff','same.diff','diff.same','same.same')
  f <- function(){
    layout(matrix(1:9, ncol=3, nrow=3, byrow=FALSE))
    par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    plot.rule <- function(rule,...){
      plot.mets <- function(ys, lb, ub, ...){
        plot(y=ys,x= 1:4, xlab='', ylab='',
             col=topo.col[c(1,3,2,4)], pch=16, cex=2,
             las=1, ylim=range.metrics,
             cex.axis=1.2, xlim=c(0,5), yaxt="n", xaxt="n",...)
        arrows(x0=1:4, y0=lb, y1=ub,
               angle=90, length=0, code=3, col=topo.col[c(1,3,2,4)])
      }
      plot.mets(ys=res.evo$metric1[res.evo$link.rule
                  == rule],
                lb=res.evo$metric1.ci.lb[res.evo$link.rule
                  == rule],
                ub=res.evo$metric1.ci.ub[res.evo$link.rule
                  == rule],
                xaxt='n', ...)
      if(rule == "matching"){
        axis(2, at=pretty(range.metrics, n=5), cex.axis=1.5, las=2)
        mtext('Edge betweenness', side=2, line=5, cex=1.5)
      }
      if(rule == "neutral"){
        legend("topright",
               legend=c('Independent evolution','No coevolution, cospeciation',
                 'Coevolution, no cospeciation',
                 'Coevolution and cospeciation'),
               col=topo.col, pch=16, bty="n", cex=1.5)
      }
      plot.mets(ys=res.evo$metric2[res.evo$link.rule
                  == rule],
                lb=res.evo$metric2.ci.lb[res.evo$link.rule
                  == rule],
                ub=res.evo$metric2.ci.ub[res.evo$link.rule
                  == rule],
                xaxt='n', ...)
      if(rule == "matching"){
        axis(2, at=pretty(range.metrics, n=5), cex.axis=1.5, las=2)
           mtext('Random walk', side=2, line=5, cex=1.5)
      }
      plot.mets(ys=res.evo$metric3[res.evo$link.rule
                  == rule],
                lb=res.evo$metric3.ci.lb[res.evo$link.rule
                  == rule],
                ub=res.evo$metric3.ci.ub[res.evo$link.rule
                  == rule],
                xaxt='n', ...)
      if(rule == "matching"){
        axis(2, at=pretty(range.metrics, n=5), cex.axis=1.5, las=2)
        mtext('Greedy', side=2, line=5, cex=1.5)
      }
      

    }
    plot.rule("matching")
    mtext('Matching', side=3, line=40, cex=1.5)
    plot.rule("barrior", yaxt="n")
    mtext('Barrier', side=3, line=40, cex=1.5)
    plot.rule("neutral", yaxt="n")
    mtext('Neutral', side=3, line=40, cex=1.5)
    ##return(res.sim)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1,
             metric2, metric3,
             sep=""))), width=9, height=9)
  
}
