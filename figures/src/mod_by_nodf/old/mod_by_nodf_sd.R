## no longer combine different abundance simulation methods
library(RColorBrewer)

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
mod.by.nodf.plot <- function(simres, path, metric1,
                             metric2, type="metrics", eco.type="eco",
                             subset=FALSE, column=NA, case=NA){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }
  
  ## simres[[2]][,"mats"][simres[[2]][,"mats"] == "eco2"] <- "eco"
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
  res.eco <- res.sim[res.sim$mat == eco.type,]
  topo.col <- brewer.pal(4, 'Spectral')
  names(topo.col) <-
    c('diff.diff','same.diff','diff.same','same.same')

  f <- function(){
    layout(matrix(1:6, ncol=3, nrow=2, byrow= FALSE))
    par(oma=c(5,9,5,1), mar=c(0.5,0,0,0.5),mgp=c(2,1,0))
    
    plot.rule <- function(rule,...){
      plot.mod.nodf <- function(dats,...){
        with(dats[dats$link.rule == rule,],{
          plot(NODF, mod, xlab='', ylab='',
               col=topo.col[topo], pch=16, cex=2,
               las=1, ylim=range.mod, xlim=range.NODF,
               cex.axis=1.2, ...)
          arrows(x0=NODF, y0=mod.ci.lb, y1=mod.ci.ub,
                 angle=90, length=0, code=3, col=topo.col[topo])
          arrows(x0=nodf.ci.lb, x1=nodf.ci.ub, y0=mod, angle=90,
                 length=0, code=3, col=topo.col[topo])
          points(NODF, mod,
                 col=topo.col[topo], pch=16, cex=2)
        })
      }
      plot.mod.nodf(res.evo, xaxt='n',...)
      if(rule == "neutral"){
        legend("topright",
               legend=c('Independent evolution',
                 'No coevolution, cospeciation',
                 'Coevolution, no cospeciation',
                 'Coevolution and cospeciation'),
               col=topo.col, pch=16, bty="n", cex=1.3)
      }
      if(rule == "matching"){
        axis(2, at=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4),
             cex.axis=1.5, las=2)
        ## axis(2, at=pretty(range.mod, n=5, min.n=3), cex.axis=1.5, las=2)
        mtext('Without abundances', side=2, line=7, cex=1.5)
      }
      plot.mod.nodf(res.eco,  xaxt='n', ...)
      if(rule == "matching"){
        axis(2, at=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4),
             cex.axis=1.5, las=2)
        ## axis(2, at=pretty(range.mod, n=5, min.n=3), cex.axis=1.5, las=2)
      }
      axis(1, at=pretty(range.NODF, n=5, min.n=3), cex.axis=1.5)
    }
    plot.rule("matching", yaxt="n")
    mtext('Matching', side=3, line=20, cex=1.5)
    mtext('With abundances', side=2, line=7, cex=1.5)
    if(type=="metrics"){
      mtext('Modularity', side=2, line=4, adj=1.5, cex=1.5)
    } else{
      mtext('Phylogenetic interaction signal (Resource species)',
            side=2, line=3, adj=-0.5)
    }
    plot.rule("barrior", yaxt="n")
    mtext('Barrier', side=3, line=20, cex=1.5)
    if(type=="metrics"){
      mtext('Nestedness', side=1, line=4, cex=1.5)
    } else{
      mtext('Phylogenetic interaction signal (Resource seeking species)',
            side=1, line=3)
    }
    plot.rule("neutral", yaxt="n")
    mtext('Neutral', side=3, line=20, cex=1.5)
    ##  return(res.sim)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf", paste(metric1,
             metric2, eco.type,
             sep=""))), width=9, height=6)
  
}
