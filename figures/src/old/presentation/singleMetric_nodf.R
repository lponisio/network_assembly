library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, path, metric1,
                             metric2, type="metrics"){

  simres[[2]][,"mats"][simres[[2]][,"mats"] == "eco2"] <- "evo"
  res.sim <- aggregate(list(NODF= simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) mean(x, na.rm=TRUE))
  sd.nodf <- aggregate(list(CI=simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"]), 
                       function(x) sd(x, na.rm=TRUE)
                       )
  
  res.sim$nodf.ci.lb <-   res.sim$NODF + sd.nodf$CI
  res.sim$nodf.ci.ub <-   res.sim$NODF - sd.nodf$CI
  means.mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                         by= list(topo= simres[[2]][,"topo"],
                           link.rule= simres[[2]][,"link.rule"],
                           mat= simres[[2]][,"mats"]), 
                         function(x) mean(x, na.rm=TRUE))
  sd.mod <- aggregate(list(CI= simres[[1]][, metric2]), 
                      by= list(topo= simres[[2]][,"topo"],
                        link.rule= simres[[2]][,"link.rule"],
                        mat= simres[[2]][,"mats"]), 
                      function(x) sd(x, na.rm=TRUE)
                      )
  
  res.sim$mod <- means.mod$mod
  res.sim$mod.ci.lb <-   res.sim$mod + sd.mod$CI 
  res.sim$mod.ci.ub <-   res.sim$mod - sd.mod$CI
  range.NODF <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub))
  range.mod <- range(c(res.sim$mod.ci.lb, res.sim$mod.ci.ub))
  
  res.evo <- res.sim[res.sim$mat == 'evo',]
  res.eco <- res.sim[res.sim$mat == 'eco',]
  topo.col <- brewer.pal(4, 'Spectral')
  names(topo.col) <-
    c('diff.diff','same.diff','diff.same','same.same')

  f <- function(){
    layout(matrix(1:3, ncol=3, nrow=1, byrow= FALSE))
    par(oma=c(5,9,5,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0),  fg="white",
        col.axis="white")
    
    plot.rule <- function(rule,...){
      plot.mod.nodf <- function(dats,...){
        with(dats[dats$link.rule == rule,],{
          plot(y=NODF, x=1:4, xlab='', ylab='',
               col=topo.col[topo], pch=16, cex=2,
               las=1, ylim=range.NODF,
               cex.axis=1.2, xlim=c(0,5),...)
          arrows(x0=1:4, y0=nodf.ci.lb, y1=nodf.ci.ub,
                 angle=90, length=0, code=3, col=topo.col[topo])
          points(y=NODF, x=1:4,
                 col=topo.col[topo], pch=16, cex=2)
        })
      }
      plot.mod.nodf(res.evo, xaxt='n',...)
      if(rule == "matching"){
        axis(2, at=pretty(range.NODF, n=5, min.n=3), cex.axis=1.5)
      }
    }
    plot.rule("matching", yaxt="n")
    plot.rule("barrior", yaxt="n")
    plot.rule("neutral", yaxt="n")
    mtext('Neutral', side=3, line=20, cex=1.5)
  ##  return(res.sim)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1, "sd",
             sep=""))), width=9, height=4)
  
}
