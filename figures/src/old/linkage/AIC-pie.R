library(RColorBrewer)
library(plotrix)

aic.plot <- function(simres, type){

  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  f <- function(){

    res.evo <- simres[simres$mats == 'evo',]
    res.eco <- simres[simres$mats == 'eco',]
    res.eco2 <- simres[simres$mats == 'eco2',]
    topo.col <- brewer.pal(4, 'Spectral')
    names(topo.col) <-
      c('diff.diff','same.diff','diff.same','same.same')
    layout(matrix(1:3, ncol=3, nrow=1, byrow= TRUE))
    with(res.evo,{
      par(oma=c(5,7,4,1), mar=c(0.5,0,0,2),mgp=c(2,1,0))
      plot(1:2, 1:2, ylim=c(0,16), xlab='', ylab='',
           col="white", pch=16, cex=2, xaxt="n",
           las=2, xlim=c(0,3), yaxt="n")
      ##matching
      floating.pie(0.75, 2, x=as.numeric(res.evo[5,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick", "firebrick4"))
      floating.pie(0.75, 6, x=as.numeric(res.evo[6,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(0.75, 10, x=as.numeric(res.evo[7,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(0.75, 14, x=as.numeric(res.evo[8,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))

      ##barrier

      floating.pie(2.25, 2, x=as.numeric(res.evo[1,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick3",
                   "firebrick4"), edges=)
      floating.pie(2.25, 6, x=as.numeric(res.evo[2,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(2.25, 10, x=as.numeric(res.evo[3,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(2.25, 14, x=as.numeric(res.evo[4,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))
      mtext('Evolutionary', side=3, line=2.5)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
      
      abline(v=1.5, lty="dashed")
      
    })
    with(res.eco,{
      plot(1:2, 1:2, ylim=c(0,16), xlab='', ylab='',
           col="white", pch=16, cex=2, xaxt="n",
           las=2, xlim=c(0,3), yaxt="n")
      floating.pie(0.75, 2, x=as.numeric(res.eco[5,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick3", "firebrick4"))
      floating.pie(0.75, 6, x=as.numeric(res.eco[6,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(0.75, 10, x=as.numeric(res.eco[7,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(0.75, 14, x=as.numeric(res.eco[8,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))

      ##barrier

      floating.pie(2.25, 2, x=as.numeric(res.eco[1,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick3", "firebrick4"))
      floating.pie(2.25, 6, x=as.numeric(res.eco[2,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(2.25, 10, x=as.numeric(res.eco[3,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(2.25, 14, x=as.numeric(res.eco[4,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))
      
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Random abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
                                       
      abline(v=1.5, lty="dashed")
                                      
    })
    with(res.eco2,{
      plot(1:2, 1:2, ylim=c(0,16), xlab='', ylab='',
           col="white", pch=16, cex=2, xaxt="n",
           las=2, xlim=c(0,3), yaxt="n")
      floating.pie(0.75, 2, x=as.numeric(res.eco2[5,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick3", "firebrick4"))
      floating.pie(0.75, 6, x=as.numeric(res.eco2[6,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(0.75, 10, x=as.numeric(res.eco2[7,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(0.75, 14, x=as.numeric(res.eco2[8,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))

      ##barrier

      floating.pie(2.25, 2, x=as.numeric(res.eco2[1,4:6]),
                   radius=0.5, col=c("firebrick1", "firebrick3", "firebrick4"))
      floating.pie(2.25, 6, x=as.numeric(res.eco2[2,4:6]),
                   radius=0.5, col=c("palegreen1", "palegreen3", "palegreen4"))
      floating.pie(2.25, 10, x=as.numeric(res.eco2[3,4:6]), radius =
                   0.5, col=c("darkorange1", "darkorange3", "darkorange4"))
      floating.pie(2.25, 14, x=as.numeric(res.eco2[4,4:6]), radius=0.5,
                   col=c("dodgerblue1", "dodgerblue3", "dodgerblue4"))
      mtext('Evo-ecological', side=3, line=2.5)
      mtext('Phylo abundances', side=3, line=1)
      mtext('Matching', side=1, line=2, adj=0.25)
      mtext('Barrier', side=1, line=2, adj=0.75)
                                       
      abline(v=1.5, lty="dashed")
                                        
    })
  }
  path <- "~/Dropbox/network_assembly/figures/coal/SingleMetrics"

  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste("aic",type,
             sep="_"))), width=9, height=5)

}

