library(RColorBrewer)
## no longer have separte simres and sim.res.samp objects
progress.plot.cor <- function(simres,  met1, met2, path, at=0.7){

  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  f <- function(){
    res.sim <- aggregate(list(metric1= simres[[1]][, met1],
                              metric2 = simres[[1]][, met2]),
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"],
                           abund= simres[[2]][,"abund"]), 
                         function(x) mean(x, na.rm=TRUE))
    res.sim$means <- mapply(function(a,b) mean(a,b, na.rm=TRUE),
                              a=res.sim[, "metric1"],
                              b=res.sim[, "metric2"])
    sd.metric <- aggregate(list(sd1= simres[[1]][, met1],
                                sd2 = simres[[1]][, met2]),
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"],
                             abund= simres[[2]][,"abund"]), 
                           function(x) sd(x, na.rm=TRUE))
    res.sim$sds <- mapply(function(a,b) mean(a,b, na.rm=TRUE),
                            a=sd.metric[,"sd1"],
                            b=sd.metric[,"sd2"])
    res.sim$sd.lb <-   res.sim$means + res.sim$sds
    res.sim$sd.ub <-   res.sim$means - res.sim$sds
    res.sim$eps <- rep(c(0.3,0.1,0.2,0), nrow(res.sim)/4)
    ranks <- cbind(unique(res.sim$mat), c(1, 2, 3, 4, 5, 6, 2))
    res.sim$case <- as.numeric(ranks[,2][match(res.sim$mat, ranks[,1])])
    res.sim <- res.sim[order(res.sim$case),] 
    range.metric <- range(c(res.sim$sd.lb, res.sim$sd.ub),
                          na.rm=TRUE)
    res.eco2 <- res.sim[res.sim$mat == "evo" |
                        res.sim$abund == "abund2",]
    
    eco2.bar <- res.eco2[res.eco2$link.rule =="barrior",]
    eco2.match <- res.eco2[res.eco2$link.rule =="matching",]
    eco2.neu <- res.eco2[res.eco2$link.rule =="neutral",]
    
    res.eco <- res.sim[res.sim$mat == "evo" |
                       res.sim$abund == "abund1",]
    
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
             las=2, type="o", pch=16, cex=1.5, lwd=1.5,
             xlim=c(2,6.4), ...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case + eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        mtext(rule, side=3, line=2, cex=1.5)
        ## mtext("Abundace", side=3, line=1)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=jitter(means), ylim=range.metric, xlab='', 
             col=topo.col[topo], xaxt="n",
             las=2, type="o",pch=16, cex=1.5, lwd=1.5,  xlim=c(2,6.4),
             ...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case + eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        mtext('100', side=1, line=2, adj=0.05, cex=1.2)
        ## mtext('80', side=1, line=2, adj=0.325, cex=0.8)
        mtext('60', side=1, line=2, adj=0.55, cex=1.2)
       ##  mtext('40', side=1, line=2, adj=0.775, cex=0.8)
        mtext('20', side=1, line=2, adj=1, cex=1.2)
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
    plotevo1 <- function(dats,...){
      with(dats,{
        plot(x= case +eps, y=jitter(means),
             col=topo.col[topo], 
             las=2, pch=16, cex=1.5, xlab="", ylab="", xaxt="n",
             ylim= range.metric, xlim=c(0.95,1.4),...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        ##mtext("No abundance", side=3, line=1)
      })
    }
    plotevo2 <- function(dats,...){
      with(dats,{
        plot(x= case +eps, y=jitter(means),
             col=topo.col[topo], 
             las=2, pch=16, cex=1.5, xlab="", ylab="", xaxt="n",
             ylim= range.metric, xlim=c(0.95,1.4),...)
        arrows(y0=sd.lb, y1=sd.ub, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    layout(matrix(1:10, ncol=5, nrow=2, byrow= TRUE),
           widths=c(1,5,1,5,5))
    par(oma=c(6,13,4,1), mar=c(1,0,0.5,0.5),mgp=c(2,1,0), fg="black")
    plotevo1(eco.match[ eco.match$mat
                       == "evo",], yaxt="n")
    axis(2, at=pretty(range.metric, n=4, min.n=3), cex.axis=1.5, las=2)
    mtext('Unstructured', side=2, line=11, cex=1.5)
    mtext('abundances', side=2, line=9, cex=1.5)
    plotrow1(eco.match[eco.match$topo == "same.same" & eco.match$mat !=
                       "evo",], 'Matching', yaxt="n")
    points1(eco.match[eco.match$topo == "same.diff" & eco.match$mat !=
                      "evo",])
    points1(eco.match[eco.match$topo == "diff.diff" & eco.match$mat !=
                      "evo",])
    points1(eco.match[eco.match$topo == "diff.same" & eco.match$mat !=
                      "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo1(eco.bar[eco.bar$mat
                     == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow1(eco.bar[eco.bar$topo == "same.same" & eco.bar$mat !=
                     "evo",], 'Barrier', yaxt="n")
    points1(eco.bar[eco.bar$topo == "same.diff" & eco.bar$mat !=
                    "evo",])
    points1(eco.bar[eco.bar$topo == "diff.diff" & eco.bar$mat !=
                    "evo",])
    points1(eco.bar[eco.bar$topo == "diff.same" & eco.bar$mat !=
                    "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotrow1(eco.neu[eco.neu$topo == "same.same" & eco.neu$mat
                     != "evo",], 'Neutral', yaxt="n")
    points1(eco.neu[eco.neu$topo == "same.diff" & eco.neu$mat
                    != "evo",])
    points1(eco.neu[eco.neu$topo == "diff.diff" & eco.neu$mat
                    != "evo",])
    points1(eco.neu[eco.neu$topo == "diff.same" & eco.neu$mat
                    != "evo",])
    legend("topright",
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white", cex=1.7, bty="n")
    par(mar=c(1,0,0.5,0.5))
    plotevo2(eco2.match[eco2.match$mat
                        == "evo",], yaxt="n")
    axis(2, at=pretty(range.metric, n=4, min.n=3), cex.axis=1.5, las=2)
     mtext('Phylogenetically', side=2, line=11,
          cex=1.5)
    mtext('structured abundances', side=2, line=9, cex=1.5)
    mtext("Phylogenetic interaction signal", 2, line=5, at=at, cex=1.5)
    plotrow2(eco2.match[eco2.match$topo == "same.same" & eco2.match$mat !=
                        "evo",], yaxt="n")
    points1(eco2.match[eco2.match$topo == "same.diff" & eco2.match$mat !=
                       "evo",])
    points1(eco2.match[eco2.match$topo == "diff.diff" & eco2.match$mat !=
                       "evo",])
    points1(eco2.match[eco2.match$topo == "diff.same" & eco2.match$mat !=
                       "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo2(eco2.bar[eco2.bar$mat
                      == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow2(eco2.bar[eco2.bar$topo == "same.same" & eco2.bar$mat !=
                      "evo",], yaxt="n")
    mtext("Percent Sampling", 1, line=5, cex=1.5, adj=-0.4)
    points1(eco2.bar[eco2.bar$topo == "same.diff" & eco2.bar$mat !=
                     "evo",])
    points1(eco2.bar[eco2.bar$topo == "diff.diff" & eco2.bar$mat !=
                     "evo",])
    points1(eco2.bar[eco2.bar$topo == "diff.same" & eco2.bar$mat !=
                     "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotrow2(eco2.neu[eco2.neu$topo == "same.same" & eco2.neu$mat
                      != "evo",], yaxt="n")
    points1(eco2.neu[eco2.neu$topo == "same.diff" & eco2.neu$mat
                     != "evo",])
    points1(eco2.neu[eco2.neu$topo == "diff.diff" & eco2.neu$mat
                     != "evo",])
    points1(eco2.neu[eco2.neu$topo == "diff.same" & eco2.neu$mat
                     != "evo",])
  }

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("ProgCor",
             sep=""))), width=12, height=8)

}
