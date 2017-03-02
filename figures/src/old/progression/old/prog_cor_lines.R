library(RColorBrewer)

progress.plot.cor <- function(simres, simres.samp, met1, met2, path){

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
                              a=res.sim[-(129:132),"metric1"],
                              b=res.sim[-(129:132),"metric2"]),
                       rep(NA,4))
    sd.metric <- aggregate(list(sd1= simres[[1]][, met1],
                                sd2 = simres[[1]][, met2]),
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) sd(x, na.rm=TRUE))
    res.sim$sds <- c(mapply(function(a,b) mean(a,b),
                            a=sd.metric[-(129:132),"sd1"],
                            b=sd.metric[-(129:132),"sd2"]),
                     rep(NA,4))
    res.sim$sd.lb <-   res.sim$means + res.sim$sds
    res.sim$sd.ub <-   res.sim$means - res.sim$sds
    res.sim$eps <- rep(c(0.3,0.1,0.2,0), 33)
    res.sim$case <- c(rep(3,12), rep(4,12), rep(5,12), rep(6,12),
                      rep(3,12), rep(4,12), rep(5,12), rep(6,12), rep(2,24),
                      rep(1,12))
    res.sim <- res.sim[order(res.sim$case),] 
    range.metric <- range(c(res.sim$sd.lb, res.sim$sd.ub),
                          na.rm=TRUE)
    res.eco2 <- res.sim[res.sim$mat == "evo" |
                        res.sim$mat == "eco2" |
                        res.sim$mat == "abund2.samp.2" |
                        res.sim$mat   == "abund2.samp.4" |
                        res.sim$mat == "abund2.samp.6" |
                        res.sim$mat == "abund2.samp.8", ]
    eco2.bar <- res.eco2[res.eco2$link.rule =="barrior",]
    eco2.match <- res.eco2[res.eco2$link.rule =="matching",]
    eco2.neu <- res.eco2[res.eco2$link.rule =="neutral",]
    res.eco <- res.sim[res.sim$mat == "evo" |
                       res.sim$mat == "eco" |
                       res.sim$mat == "abund.samp.2" |
                       res.sim$mat   == "abund.samp.4" |
                       res.sim$mat == "abund.samp.6" |
                       res.sim$mat == "abund.samp.8", ]
    eco.bar <- res.eco[res.eco$link.rule =="barrior",]
    eco.match <- res.eco[res.eco$link.rule =="matching",]
    eco.neu <- res.eco[res.eco$link.rule =="neutral",]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    plotrow1 <- function(dats, rule, ...){
      with(dats,{
        plot(x= case + eps, y=jitter(means), ylim=range.metric, xlab='', 
             col=topo.col[topo], xaxt="n",
             las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case + eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        mtext(rule, side=3, line=2)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=jitter(means), ylim=range.metric, xlab='', 
             col=topo.col[topo], xaxt="n",
             las=2, type="o",pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case + eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        mtext('Evo', side=1, line=2, adj=0)
        mtext('100%', side=1, line=2, adj=0.25)
        mtext('60%', side=1, line=2, adj=0.65)
        mtext('20%', side=1, line=2, adj=1)
      })
    }
    points1 <- function(dats,...){
      with(dats,{
        points(x= case +eps, y=jitter(means),
               col=topo.col[topo], 
               las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    layout(matrix(1:6, ncol=3, nrow=2, byrow= TRUE))
    par(oma=c(5,7,4,1), mar=c(1,1,0.5,0.5),mgp=c(2,1,0))
    plotrow1(eco.match[eco.match$topo == "same.same",], 'Matching')
    mtext("Random abundances", 2, line=3)
    points1(eco.match[eco.match$topo == "same.diff",])
    points1(eco.match[eco.match$topo == "diff.diff",])
    points1(eco.match[eco.match$topo == "diff.same",])
    plotrow1(eco.bar[eco.match$topo == "same.same",], 'Barrier', yaxt="n")
    points1(eco.bar[eco.match$topo == "same.diff",])
    points1(eco.bar[eco.match$topo == "diff.diff",])
    points1(eco.bar[eco.match$topo == "diff.same",])
    plotrow1(eco.neu[eco.match$topo == "same.same",], 'Neutral', yaxt="n")
    points1(eco.neu[eco.match$topo == "same.diff",])
    points1(eco.neu[eco.match$topo == "diff.diff",])
    points1(eco.neu[eco.match$topo == "diff.same",])
    legend("topright",
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white", cex=1.1)
    plotrow2(eco2.match[eco2.match$topo == "same.same",])
    points1(eco2.match[eco2.match$topo == "same.diff",])
    points1(eco2.match[eco2.match$topo == "diff.diff",])
    points1(eco2.match[eco2.match$topo == "diff.same",])
    mtext("Phylogenetic interaction signal", 2, line=5, adj=-30)
    mtext("Phylo abundances", 2, line=3)
    plotrow2(eco2.bar[eco2.match$topo == "same.same",], yaxt="n")
    points1(eco2.bar[eco2.match$topo == "same.diff",])
    points1(eco2.bar[eco2.match$topo == "diff.diff",])
    points1(eco2.bar[eco2.match$topo == "diff.same",])
    plotrow2(eco2.neu[eco2.match$topo == "same.same",], yaxt="n")
    points1(eco2.neu[eco2.match$topo == "same.diff",])
    points1(eco2.neu[eco2.match$topo == "diff.diff",])
    points1(eco2.neu[eco2.match$topo == "diff.same",])
  }

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("lines","cor",
             sep=""))), width=9, height=6)

}
