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
                      rep(3,12), rep(4,12), rep(5,12), rep(6,12), rep(2,24),
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
        backcols <- hsv(c(0.2,0.4,0.6), alpha=0.1)
        plot(x= case + eps, y=jitter(metric1), ylim=range.metric1, xlab='', 
             col="white", xaxt="n",
             las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        rect(xleft= 0, xright= 1.7, ybottom=-100, ytop=160,
             col=backcols[1], border=NA)
        rect(xleft= 1.7, xright= 2.7, ybottom=-100, ytop=160,
             col=backcols[2], border=NA)
        rect(xleft= 2.7, xright= 6.5, ybottom=-100, ytop=160,
             col=backcols[3], border=NA)
        points(x= case +eps, y=jitter(metric1),
               col=topo.col[topo], 
               las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case + eps, angle=90,
               length=0, code=0, col=topo.col[topo])
        mtext(rule, side=3, line=2)
      })
    }
    plotrow2 <- function(dats, ...){
      with(dats,{
        backcols <- hsv(c(0.2,0.4,0.6), alpha=0.1)
        plot(x= case + eps, y=jitter(metric2), ylim=range.metric2, xlab='', 
             col="white", xaxt="n",
             las=2, type="o",pch=16, cex=1.5, lwd=1.5,...)
        rect(xleft= 0, xright= 1.7, ybottom=-100, ytop=160,
             col=backcols[1], border=NA)
        rect(xleft= 1.7, xright= 2.7, ybottom=-100, ytop=160,
             col=backcols[2], border=NA)
        rect(xleft= 2.7, xright= 6.5, ybottom=-100, ytop=160,
             col=backcols[3], border=NA)
        points(x= case +eps, y=jitter(metric2),
               col=topo.col[topo], 
               las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case + eps, angle=90,
               length=0, code=0,  col=topo.col[topo])
        mtext('Evo', side=1, line=2, adj=0)
        mtext('100%', side=1, line=2, adj=0.25)
        mtext('60%', side=1, line=2, adj=0.65)
        mtext('20%', side=1, line=2, adj=1)
      })
    }
    points1 <- function(dats,...){
      with(dats,{
        points(x= case +eps, y=jitter(metric1),
               col=topo.col[topo], 
               las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb1, y1=sd.ub1, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    points2 <- function(dats,...){
      with(dats,{
        points(x= case +eps, y=jitter(metric2),
               col=topo.col[topo], 
               las=2, type="o", pch=16, cex=1.5, lwd=1.5,...)
        arrows(y0=sd.lb2, y1=sd.ub2, x0=case +eps, angle=90,
               length=0, code=0, col=topo.col[topo])
      })
    }
    layout(matrix(1:6, ncol=3, nrow=2, byrow= TRUE))
    par(oma=c(5,7,4,1), mar=c(1,1,0.5,0.5),mgp=c(2,1,0))
    plotrow1(res.match[res.match$topo == "same.same",], 'Matching')
    mtext("Nestedness", 2, line=3)
    points1(res.match[res.match$topo == "same.diff",])
    points1(res.match[res.match$topo == "diff.diff",])
    points1(res.match[res.match$topo == "diff.same",])
    plotrow1(res.bar[res.match$topo == "same.same",], 'Barrier', yaxt="n")
    points1(res.bar[res.match$topo == "same.diff",])
    points1(res.bar[res.match$topo == "diff.diff",])
    points1(res.bar[res.match$topo == "diff.same",])
    plotrow1(res.neu[res.match$topo == "same.same",], 'Neutral', yaxt="n")
    points1(res.neu[res.match$topo == "same.diff",])
    points1(res.neu[res.match$topo == "diff.diff",])
    points1(res.neu[res.match$topo == "diff.same",])
    legend("topright",
           legend=c('Neutral evolution','No coevolution, cospeciation',
             'Coevolution, no cospeciation','Coevolution and cospeciation'),
           col=topo.col, pch=16, bg="white", cex=1.1)
    plotrow2(res.match[res.match$topo == "same.same",])
    points2(res.match[res.match$topo == "same.diff",])
    points2(res.match[res.match$topo == "diff.diff",])
    points2(res.match[res.match$topo == "diff.same",])
    mtext("Modularity", 2, line=3)
    plotrow2(res.bar[res.match$topo == "same.same",], yaxt="n")
    points2(res.bar[res.match$topo == "same.diff",])
    points2(res.bar[res.match$topo == "diff.diff",])
    points2(res.bar[res.match$topo == "diff.same",])
    plotrow2(res.neu[res.match$topo == "same.same",], yaxt="n")
    points2(res.neu[res.match$topo == "same.diff",])
    points2(res.neu[res.match$topo == "diff.diff",])
    points2(res.neu[res.match$topo == "diff.same",])
    
  }
  
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("lines",met1, met2,
             sep=""))), width=9, height=6)

}
