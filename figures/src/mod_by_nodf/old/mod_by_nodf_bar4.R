library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mets4.bar <- function(simres, path, metric1,
                      metric2, metric3, metric4, 
                      subset=FALSE, column=NA, case=NA){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }

  ## metric 1
  res.sim <- aggregate(list(met1= simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) mean(x, na.rm=TRUE))
  sd.met1 <- aggregate(list(CI=simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) sd(x, na.rm=TRUE))
  res.sim$sd.met1 <- sd.met1$CI
  res.sim$met1.ci.lb <-   res.sim$met1 + sd.met1$CI
  res.sim$met1.ci.ub <-   res.sim$met1 - sd.met1$CI
  range.met1 <- range(c(res.sim$met1.ci.lb, res.sim$met1.ci.ub),
                      na.rm=TRUE)

  ## metric 2
  means.met2 <- aggregate(list(met2= simres[[1]][, metric2]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) mean(x, na.rm=TRUE))
  sd.met2 <- aggregate(list(CI= simres[[1]][, metric2]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) sd(x, na.rm=TRUE))
  res.sim$met2 <- means.met2$met2
  res.sim$sd.met2 <- sd.met2$CI
  res.sim$met2.ci.lb <-   res.sim$met2 + sd.met2$CI 
  res.sim$met2.ci.ub <-   res.sim$met2 - sd.met2$CI
  range.met2 <- range(c(res.sim$met2.ci.lb, res.sim$met2.ci.ub),
                      na.rm=TRUE)

  ## metric 3
  means.met3 <- aggregate(list(met3= simres[[1]][, metric3]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) mean(x, na.rm=TRUE))
  sd.met3 <- aggregate(list(CI= simres[[1]][, metric3]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) sd(x, na.rm=TRUE))
  res.sim$met3 <- means.met3$met3
  res.sim$sd.met3 <- sd.met3$CI
  res.sim$met3.ci.lb <-   res.sim$met3 + sd.met3$CI 
  res.sim$met3.ci.ub <-   res.sim$met3 - sd.met3$CI
  range.met3 <- range(c(res.sim$met3.ci.lb, res.sim$met3.ci.ub),
                      na.rm=TRUE)

  ## metric 4
  means.met4 <- aggregate(list(met4= simres[[1]][, metric4]), 
                          by= list(topo= simres[[2]][,"topo"],
                            link.rule= simres[[2]][,"link.rule"],
                            mat= simres[[2]][,"mats"]), 
                          function(x) mean(x, na.rm=TRUE))
  sd.met4 <- aggregate(list(CI= simres[[1]][, metric4]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) sd(x, na.rm=TRUE))
  res.sim$met4 <- means.met4$met4
  res.sim$sd.met4 <- sd.met4$CI
  res.sim$met4.ci.lb <-   res.sim$met4 + sd.met4$CI 
  res.sim$met4.ci.ub <-   res.sim$met4 - sd.met4$CI
  range.met4 <- range(c(res.sim$met4.ci.lb, res.sim$met4.ci.ub),
                      na.rm=TRUE)
  
  res.evo <- res.sim[res.sim$mat == 'evo',]
  topo.col <-  brewer.pal(11, 'Spectral')[c(1,3,10,11)]
  names(topo.col) <-
    c('diff.diff','same.diff','diff.same','same.same')

  f <- function(){
    plot.met1 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met1, xlab='', ylab='',
                      col=topo.col, pch=16, cex=2,
                      las=1, ylim=range.met1,
                      cex.axis=1.2, ...)
        abline(h=0, col="gray")
        segments(mp, met1 - sd.met1, mp, met1 + sd.met1,
                 lwd=2)
      })
    }
    plot.met2 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met2, xlab='', ylab='',
                      col=topo.col, pch=16, cex=2,
                      las=1, ylim=range.met2,
                      cex.axis=1.2, ...)
        abline(h=0, col="gray")
        segments(mp, met2 - sd.met2, mp, met2 + sd.met2,
                 lwd=2)
      })
    }
    plot.met3 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met3, xlab='', ylab='',
                      col=topo.col, pch=16, cex=2,
                      las=1, ylim=range.met3,
                      cex.axis=1.2, ...)
        abline(h=0, col="gray")
        segments(mp, met3 - sd.met3, mp, met3 + sd.met3,
                 lwd=2)
      })
    }
    plot.met4 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met4, xlab='', ylab='',
                      col=topo.col, pch=16, cex=2,
                      las=1, ylim=range.met4,
                      cex.axis=1.2,...)
        abline(h=0, col="gray")
        segments(mp, met4 - sd.met4, mp, met4 + sd.met4,
                 lwd=2)
        labels <- c('Independent \n evolution',
                    'No coevolution, \n cospeciation',
                    'Coevolution, \n no cospeciation',
                    'Coevolution \n and cospeciation')
        text(mp, par("usr")[3] + 0.1,
             srt = 45, adj = 1,
             labels = labels,
             xpd = NA,
             cex=1)
      })
    }
    ## layout(matrix(1:12, ncol=3, byrow= FALSE))
    ## par(oma=c(5,9,5,1), mar=c(0.5,0,2,0.5),mgp=c(2,1,0))
    ## layout(matrix(1:8, ncol=2, byrow= FALSE))
    ## par(oma=c(0,3,0.5,0),
    ##     mar=c(5,2,0,0.5),
    ##     mgp=c(2,1,0))
   layout(matrix(1:8, 4, 2, byrow=FALSE),
           widths=rep(1,8), heights=rep(1,8))
    par(oma=c(7,3,2,0.1), mar=c(0.5,3,1.5,0.1),
        mgp=c(0,0.15,0), tcl=0, cex.axis=0.8)
    
    plot.met1(res.evo, "qual", xaxt='n')
    mtext('Nestedness', side=2, line=3, cex=1)
    mtext('Qualitative', side=3, line=1, cex=1.5)
    plot.met2(res.evo, "qual", xaxt='n')
    mtext('Modularity', side=2, line=3, cex=1)
    plot.met3(res.evo, "qual", xaxt='n')
    mtext('Connectance', side=2, line=3, cex=1)
    plot.met4(res.evo, "qual", xaxt='n')
    mtext('Phylo Signal', side=2, line=3, cex=1)
    
    plot.met1(res.evo, "quan", xaxt='n', yaxt='n')
    mtext('Quantitative', side=3, line=1, cex=1.5)
    plot.met2(res.evo, "quan", xaxt='n', yaxt='n')
    plot.met3(res.evo, "quan", xaxt='n',  yaxt='n')
    plot.met4(res.evo, "quan", xaxt='n',  yaxt='n')
  }
  
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("bar", metric1,
             metric2, metric3, metric4,
             sep=""))), height=8, width=5)

}
