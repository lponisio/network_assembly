library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, fig.height, path){
  f <- function(){
    res.sim <- aggregate(list(NODF= simres[[1]][, "NODF"]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
    sd.nodf <- aggregate(list(NODF= simres[[1]][, "NODF"]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) sd(x, na.rm=TRUE))
    res.sim$sd.NODF <- sd.nodf$NODF
    means.mod <- aggregate(list(mod= simres[[1]][, "modularity"]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) mean(x, na.rm=TRUE))
    sd.mod <- aggregate(list(mod= simres[[1]][, "modularity"]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"]), 
                        function(x) sd(x, na.rm=TRUE))
    res.sim$mod <- means.mod$mod
    res.sim$sd.mod <- sd.mod$mod
    res.evo <- res.sim[res.sim$mat == 'evo',]
    res.eco <- res.sim[res.sim$mat == 'eco',]

    layout(matrix(1:2, ncol=2, nrow=1))
    topo.pch <- c(15:18)
    names(topo.pch) <- c('diff.diff','same.diff','diff.same','same.same')
    lrule.col <- c('dodgerblue','violet','springgreen')
    names(lrule.col) <- c('matching','barrior','neutral')

    with(res.evo,{
      par(mar=c(3,4,2,1))
      plot(NODF,mod, ylim=c(-.05,0.5), xlab='Nestedness', ylab="Modularity",
           col=lrule.col[link.rule], pch=topo.pch[topo], cex=2)
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3) 
    })

    with(res.eco,{
      par(mar=c(3,2,2,2))
      plot(NODF,mod, ylim=c(0,0.5), xlab='Nestedness', ylab="", yaxt="n",
           col=lrule.col[link.rule], pch=topo.pch[topo], cex=2, xlim=c(0,100))
      arrows(x0=NODF, y0=mod-sd.mod, y1=mod+sd.mod,
             angle=90, length=0.05, code=3)
      arrows(x0=NODF-sd.NODF, x1=NODF+sd.NODF, y0=mod, angle=90,
             length=0.05, code=3) 
    })

    legend("topright", legend=c("Matching", "Barrior", "Neutral"),
           col=lrule.col, pch=16, bty="n")
  }
  path <- "~/Dropbox/network_assembly/figures/coal/"

  pdf.f(f, file= file.path(path, "nodfmod.pdf"), width=7, height=3.5)
  
}

