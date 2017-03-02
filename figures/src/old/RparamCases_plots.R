library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
Rparam.plot <- function(simres, metric, fig.height, ymax, case){
  f <- function(){
    means.sim <- aggregate(list(metrics= simres[[1]][, metric]), 
                           by= list(topo= simres[[2]][,"topo"],
                             link.rule= simres[[2]][,"link.rule"],
                             mat= simres[[2]][,"mats"]), 
                           function(x) mean(x, na.rm=TRUE))
    sd.sim <- aggregate(list(metrics= simres[[1]][, metric]), 
                        by= list(topo= simres[[2]][,"topo"],
                          link.rule= simres[[2]][,"link.rule"],
                          mat= simres[[2]][,"mats"]), 
                        function(x) sd(x, na.rm=TRUE))
    means.sim$sd <- sd.sim$metrics         

    if(case=='evo'){
      means.sim <- means.sim[means.sim$mat == 'evo',]
    }

    if(case=='eco'){
      means.sim <- means.sim[means.sim$mat == 'eco',]
    }

    layout(matrix(1, 1, 1))
    topo.ord <- 1:4
    names(topo.ord) <- c('diff.diff','same.diff','diff.same','same.same')
    lrule.col <- c('dodgerblue','violet','springgreen')
    names(lrule.col) <- c('matching','barrior','neutral')

    with(means.sim,{
      par(mar=c(2,4,2,1))
      plot(topo.ord[topo],metrics, ylim=c(0,ymax),
           xaxt='n', xlab='', ylab=metric,
           col=lrule.col[link.rule], pch=16, cex=2)
      arrows(x0=topo.ord[topo], y0=metrics-sd, y1=metrics+sd,
             angle=90, length=0.05, code=3)
    })
  }

  path <- "~/Dropbox/network_assembly/figures/Rparam"
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(case, metric,
             sep=""))), width=3.5, height=3.5)
  
}

