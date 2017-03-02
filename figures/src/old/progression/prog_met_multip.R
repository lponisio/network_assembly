library(RColorBrewer)
## no longer have separte simres and sim.res.samp objects
progress.plot <- function(simres, met1, met2, path){
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
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.metric <- aggregate(list(sd1=simres[[1]][, met1],
                                sd2=simres[[1]][, met2]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) sd(x,na.rm=TRUE))
    res.sim$sd.lb1 <-   res.sim$metric1 + sd.metric$sd1
    res.sim$sd.ub1 <-   res.sim$metric1 - sd.metric$sd1
    res.sim$sd.lb2 <-   res.sim$metric2 + sd.metric$sd2
    res.sim$sd.ub2 <-   res.sim$metric2 - sd.metric$sd2
    res.sim$eps <- rep(c(0.3,0.1,0.2,0), nrow(res.sim)/4)
    ranks <- cbind(unique(res.sim$mat), c(2, 2, 1, 3, 4, 5, 6))

    res.sim$case <- as.numeric(ranks[,2][match(res.sim$mat, ranks[,1])])
    res.sim <- res.sim[order(res.sim$case),] 
    res.sim <- res.sim[res.sim$mat != "eco",]
    res.sim <- res.sim[order(res.sim$case),]
    range.metric1 <- range(c(res.sim$sd.lb1, res.sim$sd.ub1),
                           na.rm=TRUE)
    range.metric2 <- range(c(res.sim$sd.lb2, res.sim$sd.ub2),
                           na.rm=TRUE)
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    res.bar <- res.sim[res.sim$link.rule =="barrior",]
    res.match <- res.sim[res.sim$link.rule =="matching",]
    res.neu <- res.sim[res.sim$link.rule =="neutral",]
    plotrow1 <- function(dats,rule, ...){
      with(dats,{
        plot(x= case + eps, y=jitter(metric1), ylim=range.metric1,
             xlab='', 
             col= topo.col[topo], xaxt="n",
             las=2, type="o", pch=16, cex=1.5, lwd=1.5,
             xlim=c(2,6.4),...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case + eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        mtext(rule, side=3, line=2, cex=1.5)
       ## mtext("Abundance", side=3, line=1)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        plot(x= case + eps, y=jitter(metric2), ylim=range.metric2,
             xlab='', 
             col= topo.col[topo], xaxt="n",
             las=2, type="o",pch=16, cex=1.5, lwd=1.5, xlim=c(2,6.4),...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case + eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        mtext('100', side=1, line=2, adj=0.05, cex=1.2)
       # mtext('80', side=1, line=2, adj=0.325, cex=0.8)
        mtext('60', side=1, line=2, adj=0.55, cex=1.2)
       # mtext('40', side=1, line=2, adj=0.775, cex=0.8)
        mtext('20', side=1, line=2, adj=1, cex=1.2)
      })
    }
    points1 <- function(dats,...){
      with(dats,{
        points(x= case +eps, y=jitter(metric1),
               col=topo.col[topo], 
               type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    points2 <- function(dats,...){
      with(dats,{
        points(x= case +eps, y=jitter(metric2),
               col=topo.col[topo], 
               type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    plotevo1 <- function(dats,...){
      with(dats,{
        plot(x= case +eps, y=jitter(metric1),
             col=topo.col[topo], 
             las=2, pch=16, cex=1.5, xlab="", ylab="", xaxt="n",
             ylim= range.metric1, xlim=c(0.95,1.4), yaxt="n",...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
       ## mtext("No abundance", side=3, line=1)
      })
    }
    plotevo2 <- function(dats,...){
      with(dats,{
        plot(x= case +eps, y=jitter(metric2),
             col=topo.col[topo], 
             las=2, pch=16, cex=1.5, xlab="", ylab="", xaxt="n",
             ylim= range.metric2, xlim=c(0.95,1.4), yaxt="n",...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    ##browser()
    layout(matrix(1:12, ncol=6, nrow=2, byrow= TRUE), widths=c(1,5,1,5,1,5))
    par(oma=c(6,8,4,1), mar=c(1,0,0.5,0.5),mgp=c(2,1,0), fg="black")
    plotevo1(res.match[ res.match$mat
                       == "evo",])
    mtext("Nestedness", 2, line=5, cex=1.5)
     axis(2, at=pretty(range.metric1, n=4, min.n=3), cex.axis=1.5, las=2)
    plotrow1(res.match[res.match$topo == "same.same" & res.match$mat !=
                       "evo",], 'Matching', yaxt="n")
    points1(res.match[res.match$topo == "same.diff" & res.match$mat !=
                      "evo",])
    points1(res.match[res.match$topo == "diff.diff" & res.match$mat !=
                      "evo",])
    points1(res.match[res.match$topo == "diff.same" & res.match$mat !=
                      "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo1(res.bar[res.bar$mat
                     == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow1(res.bar[res.bar$topo == "same.same" & res.bar$mat !=
                     "evo",], 'Barrier', yaxt="n")
    points1(res.bar[res.bar$topo == "same.diff" & res.bar$mat !=
                    "evo",])
    points1(res.bar[res.bar$topo == "diff.diff" & res.bar$mat !=
                    "evo",])
    points1(res.bar[res.bar$topo == "diff.same" & res.bar$mat !=
                    "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo1(res.neu[res.neu$mat
                     == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow1(res.neu[res.neu$topo == "same.same" & res.neu$mat
                     != "evo",], 'Neutral', yaxt="n")
    points1(res.neu[res.neu$topo == "same.diff" & res.neu$mat
                    != "evo",])
    points1(res.neu[res.neu$topo == "diff.diff" & res.neu$mat
                    != "evo",])
    points1(res.neu[res.neu$topo == "diff.same" & res.neu$mat
                    != "evo",])
    legend("topright",
           legend=c('Independent evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white", cex=1.7, bty="n")    
    plotevo2(res.match[res.match$mat
                       == "evo",])
    mtext("Modularity", 2, line=5, cex=1.5)
    axis(2, at=c(-10, -5, 0, 5), labels=c(-10, -5, 0, 5), cex.axis=1.5, las=2)
    plotrow2(res.match[res.match$topo == "same.same" & res.match$mat !=
                       "evo",], yaxt="n")
    points2(res.match[res.match$topo == "same.diff" & res.match$mat !=
                      "evo",])
    points2(res.match[res.match$topo == "diff.diff" & res.match$mat !=
                      "evo",])
    points2(res.match[res.match$topo == "diff.same" & res.match$mat !=
                      "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo2(res.bar[res.bar$mat
                     == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow2(res.bar[res.bar$topo == "same.same" & res.bar$mat !=
                     "evo",], yaxt="n")
    mtext("Percent Sampling", 1, line=5, cex=1.5)
    points2(res.bar[res.bar$topo == "same.diff" & res.bar$mat !=
                    "evo",])
    points2(res.bar[res.bar$topo == "diff.diff" & res.bar$mat !=
                    "evo",])
    points2(res.bar[res.bar$topo == "diff.same" & res.bar$mat !=
                    "evo",])
    par(mar=c(1,0.5,0.5,0.5))
    plotevo2(res.neu[res.neu$mat
                     == "evo",], yaxt="n")
    par(mar=c(1,0,0.5,0.5))
    plotrow2(res.neu[res.neu$topo == "same.same" & res.neu$mat
                     != "evo",], yaxt="n")
    points2(res.neu[res.neu$topo == "same.diff" & res.neu$mat
                    != "evo",])
    points2(res.neu[res.neu$topo == "diff.diff" & res.neu$mat
                    != "evo",])
    points2(res.neu[res.neu$topo == "diff.same" & res.neu$mat
                    != "evo",])
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("Prog", met1, met2,
             sep=""))), width=12, height=8)

}
