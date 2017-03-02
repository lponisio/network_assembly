library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, fig.height, path, metric1, metric2){
  f <- function(){
    sim.res <- as.data.frame(simres[[1]])
    sim.res <- cbind(sim.res, simres[[2]])

    colnames(sim.res)[colnames(sim.res) == metric1] <- "metric1"
    colnames(sim.res)[colnames(sim.res) == metric2] <- "metric2"

    sim.res <- sim.res[is.finite(sim.res$metric1) &
    is.finite(sim.res$metric2),]

    
    range.NODF <- quantile(sim.res[,"metric1"], probs= c(0.025,
                                                  0.975), na.rm=TRUE)
    range.mod <- quantile(sim.res[,"metric2"], probs= c(0.025,
                                                 0.975), na.rm=TRUE)

    ## range.mod <- range(sim.res[,"metric2"], na.rm=TRUE)
    ## range.NODF <- range(sim.res[,"metric1"], na.rm=TRUE)
    res.evo <- sim.res[sim.res$mats == 'evo',]
    res.eco <- sim.res[sim.res$mats == 'eco',]
    res.eco2 <- sim.res[sim.res$mats == 'eco2',]
    
    topo.col <- brewer.pal(4, 'Spectral')
    
    names(topo.col) <- c('diff.diff','same.diff','diff.same','same.same')

    ##**********************************************************
    ## matching
    ##**********************************************************
    ## evo
    layout(matrix(1:9, ncol=3, nrow=3, byrow= TRUE))
    with(res.evo[res.evo$link.rule == 'matching' &
                 res.evo$topo == 'same.same',],{
                   par(oma=c(5,7,4,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
                   plot(1, 1, xlab='', ylab='',
                        col="white", pch=16, cex=2, xaxt='n',
                        las=2, ylim=range.mod, xlim=range.NODF)
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                   mtext('Matching', side=2, line=5)
                   mtext('Evolutionary', side=3, line=2.2)
                 })
    with(res.evo[res.evo$link.rule == 'matching' &
                 res.evo$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'matching' &
                 res.evo$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'matching' &
                 res.evo$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ## eco
    with(res.eco[res.eco$link.rule == 'matching' &
                 res.eco$topo == 'same.same',],{
                   plot(1, 1, xlab='', ylab="", yaxt="n",
                        col="white", pch=16, cex=2, xaxt='n',
                        xlim=range.NODF, ylim=range.mod)
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                   mtext('Evo-Evological', side=3, line=2.5)
                   mtext('Random Abundances', side=3, line=1)
                 })
    with(res.eco[res.eco$link.rule == 'matching' &
                 res.eco$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'matching' &
                 res.eco$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'matching' &
                 res.eco$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ##eco 2
    with(res.eco2[res.eco2$link.rule == 'matching' &
                  res.eco2$topo == 'same.same',],{
                    plot(1, 1 ,xlab='', ylab="", yaxt="n",
                         col="white", pch=16, cex=2,
                         xaxt='n', ylim=range.mod, xlim=range.NODF)
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                    mtext('Evo-Evological', side=3, line=2.5)
                    mtext('Phylo Abundances', side=3, line=1)
                  })
    with(res.eco2[res.eco2$link.rule == 'matching' &
                  res.eco2$topo == 'same.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'matching' &
                  res.eco2$topo == 'diff.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'matching' &
                  res.eco2$topo == 'diff.same',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    ##**********************************************************
    ## barrier
    ##**********************************************************
    with(res.evo[res.evo$link.rule == 'barrior' &
                 res.evo$topo == 'same.same',],{
                   plot(1, 1, xlab='', ylab='',
                        col="white", pch=16, cex=2, xaxt='n',
                        las=2, ylim=range.mod, xlim=range.NODF)
                   kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                   try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo]), silent=TRUE)
                   mtext('Barrier', side=2, line=5)
                 })
    with(res.evo[res.evo$link.rule == 'barrior' &
                 res.evo$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'barrior' &
                 res.evo$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'barrior' &
                 res.evo$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ## eco
    with(res.eco[res.eco$link.rule == 'barrior' &
                 res.eco$topo == 'same.same',],{
                   plot(1, 1, xlab='', ylab="", yaxt="n",
                        col="white", pch=16, cex=2, xaxt='n',
                        xlim=range.NODF, ylim=range.mod)
                   kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                   try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo]), silent=TRUE)
                 })
    with(res.eco[res.eco$link.rule == 'matching' &
                 res.eco$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'barrior' &
                 res.eco$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'barrior' &
                 res.eco$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ##eco 2
    with(res.eco2[res.eco2$link.rule == 'barrior' &
                  res.eco2$topo == 'same.same',],{
                    plot(1, 1 ,xlab='', ylab="", yaxt="n",
                         col="white", pch=16, cex=2,
                         xaxt='n', ylim=range.mod, xlim=range.NODF)
                    kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                    try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo]), silent=TRUE)
                  })
    with(res.eco2[res.eco2$link.rule == 'barrior' &
                  res.eco2$topo == 'same.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'barrior' &
                  res.eco2$topo == 'diff.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'barrior' &
                  res.eco2$topo == 'diff.same',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    ##**********************************************************
    ## neutral
    ##**********************************************************
    with(res.evo[res.evo$link.rule == 'neutral' &
                 res.evo$topo == 'same.same',],{
                   plot(1, 1, xlab='', ylab='',
                        col="white", pch=16, cex=2,
                        las=2, ylim=range.mod, xlim=range.NODF)
                   kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                   try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo]), silent=TRUE)
                   mtext('Barrier', side=2, line=5)
                 })
    with(res.evo[res.evo$link.rule == 'neutral' &
                 res.evo$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'neutral' &
                 res.evo$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.evo[res.evo$link.rule == 'neutral' &
                 res.evo$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ## eco
    with(res.eco[res.eco$link.rule == 'neutral' &
                 res.eco$topo == 'same.same',],{
                   plot(1, 1, xlab='', ylab="", yaxt="n",
                        col="white", pch=16, cex=2,
                        xlim=range.NODF, ylim=range.mod)
                   kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                   try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo]), silent=TRUE)
                 })
    with(res.eco[res.eco$link.rule == 'neutral' &
                 res.eco$topo == 'same.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'neutral' &
                 res.eco$topo == 'diff.diff',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    with(res.eco[res.eco$link.rule == 'neutral' &
                 res.eco$topo == 'diff.same',],{
                   kernal <- kde2d(x=metric1, y=metric2)
                   contour(x=kernal$x, y=kernal$y, z=kernal$z,
                           nlevels=5, drawlabels=FALSE,add=TRUE,
                           col=topo.col[topo])
                 })
    ##eco 2
    with(res.eco2[res.eco2$link.rule == 'neutral' &
                  res.eco2$topo == 'same.same',],{
                    plot(1, 1 ,xlab='', ylab="", yaxt="n",
                         col="white", pch=16, cex=2,
                         ylim=range.mod, xlim=range.NODF)
                    kernal <- try(kde2d(x=metric1, y=metric2), silent=TRUE)
                    try(contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo]), silent=TRUE)
                  })
    with(res.eco2[res.eco2$link.rule == 'neutral' &
                  res.eco2$topo == 'same.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'neutral' &
                  res.eco2$topo == 'diff.diff',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    with(res.eco2[res.eco2$link.rule == 'neutral' &
                  res.eco2$topo == 'diff.same',],{
                    kernal <- kde2d(x=metric1, y=metric2)
                    contour(x=kernal$x, y=kernal$y, z=kernal$z,
                            nlevels=5, drawlabels=FALSE,add=TRUE,
                            col=topo.col[topo])
                  })
    mtext('Nestedness', side=1, line=3, outer=TRUE)
    mtext('Modularity', side=2, line=3, outer=TRUE)
  }

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1,
             metric2, "contours",
             sep=""))), width=7, height=7)
  
}

