library(RColorBrewer)

progress.plot <- function(simres, simres.samp, met1, met2, met3){

  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  f <- function(){
    simres[[1]] <- rbind(simres[[1]], simres.samp[[1]])
    simres[[2]] <- rbind(simres[[2]], simres.samp[[2]])
    res.sim <- aggregate(list(metric1= simres[[1]][, met1],
                              metric2 = simres[[1]][, met2],
                              metric3 = simres[[1]][, met3]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.metric <- aggregate(list(sd1=simres[[1]][, met1],
                                sd2=simres[[1]][, met2],
                                sd3=simres[[1]][, met3]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) sd(x,
                                          na.rm=TRUE))
    res.sim$sd.lb1 <-   res.sim$metric1 + sd.metric$sd1
    res.sim$sd.ub1 <-   res.sim$metric1 - sd.metric$sd1
    res.sim$sd.lb2 <-   res.sim$metric2 + sd.metric$sd2
    res.sim$sd.ub2 <-   res.sim$metric2 - sd.metric$sd2
    res.sim$sd.lb3 <-   res.sim$metric3 + sd.metric$sd3
    res.sim$sd.ub3 <-   res.sim$metric3 - sd.metric$sd3
    res.sim$case <- c(rep(3,12), rep(4,12), rep(3, 12), rep(4, 12),
                      rep(2,24), rep(1,12))
    res.sim$pchs <- c(rep(1,24), rep(16, 24),
                      rep(1,12), rep(16,12), rep(16,12))
    res.sim$eps <- rep(c(-0.4,-0.2,0,0.2), 21)

    res.sim <- res.sim[res.sim$mat != "eco",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.1",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.5",]
        
    range.metric1 <- range(c(res.sim$sd.lb1, res.sim$sd.ub1),
                           na.rm=TRUE)
    range.metric2 <- range(c(res.sim$sd.lb2, res.sim$sd.ub2),
                           na.rm=TRUE)
    range.metric3 <- range(c(res.sim$sd.lb3, res.sim$sd.ub3),
                           na.rm=TRUE)
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    res.bar <- res.sim[res.sim$link.rule =="barrior",]
    res.match <- res.sim[res.sim$link.rule =="matching",]
    res.neu <- res.sim[res.sim$link.rule =="neutral",]
    ##browser()
    plotrow1 <- function(dats,rule, ...){
      with(dats,{
        plot(x= case + eps, y=metric1, ylim=range.metric1, xlab='', 
             col=topo.col[topo], pch=pchs, cex=2, xaxt="n",
             las=2,...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        points(case +eps, metric1,
               col=topo.col[topo], pch=pchs, cex=2)
        abline(v=1.4, lty="dashed")
        abline(v=2.4, lty="dashed")
        abline(v=3.4, lty="dashed")
        mtext(rule, side=3, line=2.5)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=metric2, ylim=range.metric2, xlab='', 
             col=topo.col[topo], pch=pchs, cex=2, xaxt="n",
             las=2,...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case +eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        points(case +eps, metric2,
               col=topo.col[topo], pch=pchs, cex=2)
        abline(v=1.4, lty="dashed")
        abline(v=2.4, lty="dashed")
        abline(v=3.4, lty="dashed")
      })
    }
    plotrow3 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=metric3, ylim=range.metric3, xlab='', 
             col=topo.col[topo], pch=pchs, cex=2, xaxt="n",
             las=2,...)
        arrows(y0=sd.lb3, y1=sd.ub3, x0=case +eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        points(case +eps, metric3,
               col=topo.col[topo], pch=pchs, cex=2)
        abline(v=1.4, lty="dashed")
        abline(v=2.4, lty="dashed")
        abline(v=3.4, lty="dashed")
        mtext('Evo', side=1, line=2, adj=0.1)
        mtext('Eco', side=1, line=2, adj=0.4)
        mtext('Good', side=1, line=2, adj=0.7)
        mtext('Poor', side=1, line=2, adj=1)
      })
    }
    layout(matrix(1:9, ncol=3, nrow=3, byrow= TRUE))
    par(oma=c(5,7,4,1), mar=c(0.5,0.5,0,0.5),mgp=c(2,1,0))
    plotrow1(res.match, 'Matching')
    mtext("Nestedness", 2, line=3)
    plotrow1(res.bar, 'Barrier', yaxt="n")
    plotrow1(res.neu, 'Neutral', yaxt="n")
    legend("topright",
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white")
    plotrow2(res.match)
    mtext("Modularity", 2, line=3)
    plotrow2(res.bar, yaxt="n")
    plotrow2(res.neu, yaxt="n")
    plotrow3(res.match)
    mtext("Phylogenetic Interaction Signal", 2, line=3)
    plotrow3(res.bar, yaxt="n")
    plotrow3(res.neu, yaxt="n")
  }
  path <- "~/Dropbox/network_assembly/figures/coal/progression"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(met1,
             sep=""))), width=9, height=9)

}

