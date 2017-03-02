library(RColorBrewer)

mets3.bar <- function(simres,
                      path.fig,
                      metric1,
                      metric2,
                      metric3,
                      subset=FALSE,
                      column=NA,
                      case=NA,
                      met1lab='Relative \n  Modularity',
                      met2lab='Relative \n Nestedness',
                      met3lab='Connectance',
                      adj.lab=1){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }

  ## metric 1
  res.sim <- aggregate(list(met1= simres[[1]][, metric1]),
                       by= list(topo= simres[[2]][,'topo'],
                                link.rule= simres[[2]][,'link.rule'],
                                mat= simres[[2]][,'mats']),
                       function(x) mean(x, na.rm=TRUE))
  sd.met1 <- aggregate(list(CI=simres[[1]][, metric1]),
                       by= list(topo= simres[[2]][,'topo'],
                                link.rule= simres[[2]][,'link.rule'],
                                mat= simres[[2]][,'mats']),
                       function(x) sd(x, na.rm=TRUE))
  res.sim$sd.met1 <- sd.met1$CI
  res.sim$met1.ci.lb <-   res.sim$met1 + 1.96*sd.met1$CI
  res.sim$met1.ci.ub <-   res.sim$met1 - 1.96*sd.met1$CI
  range.met1 <- range(c(res.sim$met1.ci.lb, res.sim$met1.ci.ub),
                      na.rm=TRUE)

  ## metric 2
  means.met2 <- aggregate(list(met2= simres[[1]][, metric2]),
                          by= list(topo= simres[[2]][,'topo'],
                                   link.rule= simres[[2]][,'link.rule'],
                                   mat= simres[[2]][,'mats']),
                          function(x) mean(x, na.rm=TRUE))
  sd.met2 <- aggregate(list(CI= simres[[1]][, metric2]),
                       by= list(topo= simres[[2]][,'topo'],
                                link.rule= simres[[2]][,'link.rule'],
                                mat= simres[[2]][,'mats']),
                       function(x) sd(x, na.rm=TRUE))
  ## metric 3
  means.met3 <- aggregate(list(met3= simres[[1]][, metric3]),
                          by= list(topo= simres[[2]][,'topo'],
                                   link.rule= simres[[2]][,'link.rule'],
                                   mat= simres[[2]][,'mats']),
                          function(x) mean(x, na.rm=TRUE))
  sd.met3 <- aggregate(list(CI= simres[[1]][, metric3]),
                       by= list(topo= simres[[2]][,'topo'],
                                link.rule= simres[[2]][,'link.rule'],
                                mat= simres[[2]][,'mats']),
                       function(x) sd(x, na.rm=TRUE))

  ## make table for ms
  met2 <-  cbind(means.met2,
                 sd.met2$CI,
                 means.met2$met2 + sd.met2$CI,
                 means.met2$met2 - sd.met2$CI)
  met3 <-     cbind(means.met3,
                    sd.met3$CI,
                    means.met3$met3 + sd.met3$CI,
                    means.met3$met3 - sd.met3$CI)
  colnames(met2) <- colnames(met3) <- colnames(res.sim)
  out <- rbind(res.sim, met2, met3)
  out[, c(4:7)] <- apply(out[,c(4:7)], 2, round, 3)
  out$metric <- rep(c(metric1, metric2, metric3), each=8)
  out$link.rule[out$link.rule == "qual"] <- "Unweighted"
  out$link.rule[out$link.rule == "quan"] <- "Weighted"
  out <- out[, c(-3, -6, -7)]
  out <- out[, c("topo", "metric",
                 "link.rule", "met1",
                 "sd.met1")]
  write.table(out, file=file.path(path.fig, "saved/stats.txt"),
              sep="&",
              row.names=FALSE)


  ## apend to res.sim
  res.sim$met2 <- means.met2$met2
  res.sim$sd.met2 <- sd.met2$CI
  res.sim$met2.ci.lb <-   res.sim$met2 + sd.met2$CI
  res.sim$met2.ci.ub <-   res.sim$met2 - sd.met2$CI
  range.met2 <- range(c(res.sim$met2.ci.lb, res.sim$met2.ci.ub),
                      na.rm=TRUE)

  res.sim$met3 <- means.met3$met3
  res.sim$sd.met3 <- sd.met3$CI
  res.sim$met3.ci.lb <-   res.sim$met3 + sd.met3$CI
  res.sim$met3.ci.ub <-   res.sim$met3 - sd.met3$CI
  range.met3 <- range(c(res.sim$met3.ci.lb, res.sim$met3.ci.ub),
                      na.rm=TRUE)

  res.evo <- res.sim[res.sim$mat == 'evo',]

  topo.order <-  c('diff.diff','same.diff','diff.same','same.same')

  res.evo <- res.evo[order(c(topo.order, topo.order), res.evo$topo),]
  topo.col <- brewer.pal(11, 'Spectral')[c(1,3,10,11)]
  names(topo.col) <- topo.order

  f <- function(){
    plot.met1 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met1, xlab='', ylab='',
                      col=topo.col,
                      las=1,
                      ylim=range.met1,
                      cex.axis=1.2, ...)
        abline(h=0, col='gray')
        segments(mp, met1 - sd.met1, mp, met1 + sd.met1,
                 lwd=2)
      })
    }
    plot.met2 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met2, xlab='', ylab='',
                      col=topo.col,
                      las=1,
                      ylim= range.met2,
                      cex.axis=1.2, ...)
        abline(h=0, col='gray')
        segments(mp, met2 - sd.met2, mp, met2 + sd.met2,
                 lwd=2)
      })
    }
    plot.met3 <- function(dats, rule, ...){
      with(dats[dats$link.rule == rule,],{
        mp <- barplot(met3, xlab='', ylab='',
                      col=topo.col,
                      las=1,
                      ylim=range.met3,
                      cex.axis=1.2,...)
        abline(h=0, col='gray')
        segments(mp, met3 - sd.met3, mp, met3 + sd.met3,
                 lwd=2)
        labels <- c('No coevolution, \n no cospeciation',
                    'No coevolution, \n cospeciation',
                    'Coevolution, \n no cospeciation',
                    'Coevolution \n and cospeciation')
        text(mp, par('usr')[3] - adj.lab,
             srt = 45, adj = 1,
             labels = labels,
             xpd = NA,
             cex=1)
      })
    }

    layout(matrix(1:6, 3, 2, byrow=FALSE),
           widths=rep(1,8), heights=rep(1,8))
    par(oma=c(7,3,2,0.1), mar=c(0.5,3,1.5,0.1),
        mgp=c(0,0.15,0), tcl=0, cex.axis=0.8)

    plot.met1(res.evo, 'qual', xaxt='n')
    mtext(met1lab, side=2, line=3, cex=1)
    mtext("a)", side=3, line=0, cex=1, adj=1)
    mtext('Unweighted', side=3, line=1, cex=1.5)
    plot.met2(res.evo, 'qual', xaxt='n')
    mtext("c)", side=3, line=0, cex=1, adj=1)
    mtext(met2lab, side=2, line=3, cex=1)
    plot.met3(res.evo, 'qual', xaxt='n')
    mtext("e)", side=3, line=0, cex=1, adj=1)
    mtext(met3lab, side=2, line=3, cex=1)

    plot.met1(res.evo, 'quan', xaxt='n', yaxt='n')
    mtext('Weighted', side=3, line=1, cex=1.5)
    mtext("b)", side=3, line=0, cex=1, adj=1)
    plot.met2(res.evo, 'quan', xaxt='n', yaxt='n')
    mtext("d)", side=3, line=0, cex=1, adj=1)
    plot.met3(res.evo, 'quan', xaxt='n',  yaxt='n')
    mtext("f)", side=3, line=0, cex=1, adj=1)

  }

  pdf.f(f, file= file.path(path.fig, sprintf('%s.pdf',
                                         paste('bar', metric1,
                                               metric2, metric3,
                                               sep=''))), height=6, width=5)

}
