library(RColorBrewer)

progress.plot <- function(simres, simres.samp, met1, met2, path){

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
    sd.metric <- aggregate(list(sd1=simres[[1]][, met1],
                                sd2=simres[[1]][, met2]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) sd(x,
                                          na.rm=TRUE))
    res.sim$sd.lb1 <-   res.sim$metric1 + sd.metric$sd1
    res.sim$sd.ub1 <-   res.sim$metric1 - sd.metric$sd1
    res.sim$sd.lb2 <-   res.sim$metric2 + sd.metric$sd2
    res.sim$sd.ub2 <-   res.sim$metric2 - sd.metric$sd2
    res.sim$eps <- rep(c(0.3,0.1,0.2,0), 33)
    res.sim$case <- c(rep(3,12), rep(4,12), rep(5,12), rep(6,12),
                      rep(3,12), rep(4,12), rep(5,12), rep(6,12),
                      rep(2,24),
                      rep(1,12))
    res.sim <- res.sim[order(res.sim$case),] 
    res.sim <- res.sim[res.sim$mat != "eco",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.2",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.4",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.6",]
    res.sim <- res.sim[res.sim$mat != "abund.samp.8",]
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
        mtext(rule, side=3, line=2)
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
        mtext('100', side=1, line=2, adj=0.05, cex=0.8)
        mtext('80', side=1, line=2, adj=0.325, cex=0.8)
        mtext('60', side=1, line=2, adj=0.55, cex=0.8)
        mtext('40', side=1, line=2, adj=0.775, cex=0.8)
        mtext('20', side=1, line=2, adj=1, cex=0.8)
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
             ylim= range.metric1, xlim=c(0.95,1.4),...)
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
             ylim= range.metric2, xlim=c(0.95,1.4),...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    ##browser()
    layout(matrix(1:12, ncol=6, nrow=2, byrow= TRUE), widths=c(1,5,1,5,1,5))
    par(oma=c(5,7,4,1), mar=c(1,0,0.5,0.5),mgp=c(2,1,0), fg="black")
    plotevo1(res.match[ res.match$mat
                       == "evo",])
    mtext("Nestedness", 2, line=3)
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
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white", cex=1.1, bty="n")    
    plotevo2(res.match[res.match$mat
                       == "evo",])
    mtext("Modularity", 2, line=3)
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
    mtext("Percent Sampling", 1, line=4)
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
             paste("lines_multi",met1, met2,
             sep=""))), width=9, height=6)

}
