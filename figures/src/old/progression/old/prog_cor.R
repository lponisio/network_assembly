library(RColorBrewer)

progress.plot.cor <- function(simres, simres.samp, met1, met2){

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
                              metric2 = simres[[1]][, met2]),
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))

    res.sim$means <- c(mapply(function(a,b) mean(a,b),
                              a=res.sim[-(81:84),"metric1"],
                              b=res.sim[-(81:84),"metric2"]),
                       rep(NA,4))

    sd.metric <- aggregate(list(sd1= simres[[1]][, met1],
                                sd2 = simres[[1]][, met2]),
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) sd(x, na.rm=TRUE))

    res.sim$sds <- c(mapply(function(a,b) mean(a,b),
                            a=sd.metric[-(81:84),"sd1"],
                            b=sd.metric[-(81:84),"sd2"]),
                     rep(NA,4))
    
    res.sim$sd.lb <-   res.sim$means + res.sim$sds
    res.sim$sd.ub <-   res.sim$means - res.sim$sds

    res.sim$case <- c(rep(3,12), rep(4,12), rep(3, 12), rep(4, 12),
                      rep(2,24), rep(1,12))
    res.sim$eps <- rep(c(-0.4,-0.2,0,0.2), 21)
    res.sim <- res.sim[order(res.sim$case),]

    range.metric <- range(c(res.sim$sd.lb, res.sim$sd.ub),
                          na.rm=TRUE)

    res.eco2 <- res.sim[res.sim$mat == "evo" |
                         res.sim$mat == "eco2" | res.sim$mat ==
                        "abund2.samp.1"  | res.sim$mat ==
                        "abund2.samp.5", ]
    eco2.bar <- res.eco2[res.eco2$link.rule =="barrior",]
    eco2.match <- res.eco2[res.eco2$link.rule =="matching",]
    eco2.neu <- res.eco2[res.eco2$link.rule =="neutral",]

    res.eco <- res.sim[res.sim$mat == "evo" |
                       res.sim$mat == "eco" | res.sim$mat ==
                       "abund.samp.1"  | res.sim$mat ==
                       "abund.samp.5", ]
 
    eco.bar <- res.eco[res.eco$link.rule =="barrior",]
    eco.match <- res.eco[res.eco$link.rule =="matching",]
    eco.neu <- res.eco[res.eco$link.rule =="neutral",]

    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')

    plotrow1 <- function(dats,rule, ...){
      with(dats,{
        plot(x= case + eps, y=means, ylim=range.metric, xlab='', 
             col=topo.col[topo], pch=16, cex=2, xaxt="n",
             las=2,...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        points(case +eps, means,
               col=topo.col[topo], pch=16, cex=2)
        abline(v=1.4, lty="dashed")
        abline(v=2.4, lty="dashed")
        abline(v=3.4, lty="dashed")
        mtext(rule, side=3, line=2.5)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=means, ylim=range.metric, xlab='', 
             col=topo.col[topo], pch=16, cex=2, xaxt="n",
             las=2,...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case +eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        points(case +eps, means,
               col=topo.col[topo], pch=16, cex=2)
        abline(v=1.4, lty="dashed")
        abline(v=2.4, lty="dashed")
        abline(v=3.4, lty="dashed")
        mtext('Evo', side=1, line=2, adj=0.1)
        mtext('Eco', side=1, line=2, adj=0.4)
        mtext('Good', side=1, line=2, adj=0.7)
        mtext('Poor', side=1, line=2, adj=1)
      })
    }
    layout(matrix(1:6, ncol=3, nrow=2, byrow= TRUE))
    par(oma=c(5,7,4,1), mar=c(0.5,0.5,0,0.5),mgp=c(2,1,0))
    plotrow1(eco.match, 'Matching')
     mtext("Random abundances", 2, line=3)
    plotrow1(eco.bar, 'Barrier', yaxt="n")
    plotrow1(eco.neu, 'Neutral', yaxt="n")
    legend("topright",
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white")
    plotrow2(eco2.match)
    mtext("Phylogenetic interaction signal", 2, line=5, adj=12)
    mtext("Phylo abundances", 2, line=3)
    plotrow2(eco2.bar, yaxt="n")
    plotrow2(eco2.neu, yaxt="n")
  }
  path <- "~/Dropbox/network_assembly/figures/coal/progression"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("cor",
             sep=""))), width=9, height=6)

}

