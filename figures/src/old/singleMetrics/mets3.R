library(RColorBrewer)

mets3.plot <- function(simres, metric1,
                       metric2, metric3, xlabs, path.dir,
                       subset=FALSE, column=NA, case=NA){
  if(subset){
    simres[[2]] <- simres[[2]][simres[[1]][, column] == case,]
    simres[[1]] <- simres[[1]][simres[[1]][, column] == case,]
  }
  res.sim <- aggregate(list(nodf= simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"],
                         nsp = simres[[1]][, metric3]), 
                       function(x) mean(x, na.rm=TRUE))
  sd.nodf <- aggregate(list(CI=simres[[1]][, metric1]), 
                       by= list(topo= simres[[2]][,"topo"],
                         link.rule= simres[[2]][,"link.rule"],
                         mat= simres[[2]][,"mats"],
                         nsp = simres[[1]][, metric3]), 
                       function(x) sd(x, na.rm=TRUE))
  res.sim$nodf.ci.lb <-   res.sim$nodf + sd.nodf$CI
  res.sim$nodf.ci.ub <-   res.sim$nodf - sd.nodf$CI
  mod <- aggregate(list(mod= simres[[1]][, metric2]), 
                   by= list(topo= simres[[2]][,"topo"],
                     link.rule= simres[[2]][,"link.rule"],
                     mat= simres[[2]][,"mats"],
                     nsp = simres[[1]][, metric3]), 
                   function(x) mean(x, na.rm=TRUE))
  sd.mod <- aggregate(list(CI=simres[[1]][, metric2]), 
                      by= list(topo= simres[[2]][,"topo"],
                        link.rule= simres[[2]][,"link.rule"],
                        mat= simres[[2]][,"mats"],
                        nsp = simres[[1]][, metric3]), 
                      function(x) sd(x, na.rm=TRUE))
  res.sim$mod <- mod$mod
  res.sim$mod.ci.lb <-   res.sim$mod + sd.mod$CI
  res.sim$mod.ci.ub <-   res.sim$mod - sd.mod$CI

  range.nodf <- range(c(res.sim$nodf.ci.lb, res.sim$nodf.ci.ub),
                      na.rm=TRUE)
  range.mod <- range(c(res.sim$mod.ci.lb, res.sim$mod.ci.ub),
                     na.rm=TRUE)
  range.nsp <- range(res.sim$nsp, na.rm =TRUE)
  
  res.evo <- res.sim[res.sim$mat == 'evo',]
  topo.col <-  brewer.pal(11, 'Spectral')[c(1,3,10,11)]
  fill.col <- add.alpha(topo.col, alpha=0.2)
  names(topo.col) <- names(fill.col) <- 
    c('diff.diff','same.diff','diff.same','same.same')
  f <- function(){
    layout(matrix(1:4, ncol=2, nrow=2, byrow= TRUE))
    par(oma=c(5, 9, 5, 1), mar=c(0.5, 0, 0.5, 1), mgp=c(2, 1, 0))
    plot.nodf <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        not.na <- !is.na(nodf.ci.ub)
        if(topos == "same.same"){
          plot(smooth.spline(nodf[not.na] ~ nsp[not.na], df=5),
               xlab='',
               ylab='',
               col=topo.col[topo],
               pch=16,
               cex=2,
               las=1,
               ylim=range.nodf,
               xlim=range.nsp,
               yaxt="n",
               xaxt="n",
               cex.axis=1.2,
               type="l",
               lwd=2, ...)
        } else{
          points(smooth.spline(nodf[not.na] ~ nsp[not.na], df=5),
                 col=topo.col[topo], pch=16, cex=2,
                 type="l", lwd=2, ...)
        }
        out1 <- smooth.spline(nsp[not.na], nodf.ci.lb[not.na], df=5)
        out2 <- smooth.spline(nsp[not.na], nodf.ci.ub[not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y

        polygon(x=c(nsp[not.na], rev(nsp[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)
        
        ## polygon(x=c(nsp[not.na], rev(nsp[not.na])),
        ##         y=c(nodf.ci.ub[not.na], rev(nodf.ci.lb[not.na])),
        ##         col=fill.col[topo], border = NA)
        if(rule == "neutral" & topos == "diff.diff"){
          if(xlabs == "Connectance"){
            points(nodf ~ nsp,
                   col=topo.col[topo], pch="-", cex=2)
          }
          legend("topright",
                 legend=c('Independent evolution',
                   'No coevolution, cospeciation',
                   'Coevolution, no cospeciation',
                   'Coevolution and cospeciation'),
                 col=topo.col, pch=16, bty="n", cex=1.3)
        }
        if(rule == "qual" & topos == "same.same"){
          axis(2, at=pretty(range.nodf, n=5, min.n=3),
               cex.axis=1.5, las=2)
          mtext('Nestedness', side=2, line=5, cex=1.5)
        }
      })
    }
    plot.mod <- function(dats, rule, topos, ...){
      with(dats[dats$link.rule == rule & dats$topo == topos,],{
        not.na <- !is.na(mod.ci.ub)
        if(topos == "same.same"){
          plot(smooth.spline(mod ~ nsp, df=5),
               xlab='',
               ylab='',
               col=topo.col[topo],
               pch=16,
               cex=2,
               las=1,
               ylim=range.mod,
               xlim=range.nsp,
               yaxt="n",
               xaxt="n",
               cex.axis=1.2,
               type="l",
               lwd=2, ...)
          axis(1, at=pretty(range.nsp, n=5, min.n=3),
               cex.axis=1.5)
        } else{
          points(smooth.spline(mod[not.na] ~ nsp[not.na], df=5),
                 col=topo.col[topo], pch=16, cex=2,
                 type="l", lwd=2, ...)
        }
        out1 <- smooth.spline(nsp[not.na], mod.ci.lb[not.na], df=5)
        out2 <- smooth.spline(nsp[not.na], mod.ci.ub[not.na], df=5)
        y1.min <- out1$y
        y1.max <- out2$y

        polygon(x=c(nsp[not.na], rev(nsp[not.na])),
                y=c(y1.min, rev(y1.max)),
                col=fill.col[topo], border = NA)
        
        ## polygon(x=c(nsp[not.na], rev(nsp[not.na])),
        ##         y=c(mod.ci.ub[not.na], rev(mod.ci.lb[not.na])),
        ##         col=fill.col[topo], border = NA)
        if(rule == "qual" & topos == "same.same"){
          ## axis(2, at=c(-8, -4, 0, 4), labels=c(-8, -4, 0, 4),
          ##      cex.axis=1.5, las=2)
          axis(2, at=pretty(range.mod, n=5, min.n=3),
               cex.axis=1.5, las=2)
          mtext('Modularity', side=2, line=5, cex=1.5)
        }
        if(rule == "neutral" & xlabs == "Connectance"){
          points(mod ~ nsp,
                 col=topo.col[topo], pch="-", cex=2)
        }
        
      })
    }
    plot.nodf(res.evo, "qual", "same.same")
    plot.nodf(res.evo, "qual", "same.diff")
    plot.nodf(res.evo, "qual", "diff.same")
    plot.nodf(res.evo, "qual", "diff.diff")
    mtext('Qualitative', side=3, line=2, cex=1.5)
    plot.nodf(res.evo, "quan", "same.same")
    plot.nodf(res.evo, "quan", "same.diff")
    plot.nodf(res.evo, "quan", "diff.same")
    plot.nodf(res.evo, "quan", "diff.diff")
    mtext('Quantitative', side=3, line=2, cex=1.5)
    ## plot.nodf(res.evo, "neutral", "same.same")
    ## plot.nodf(res.evo, "neutral", "same.diff")
    ## plot.nodf(res.evo, "neutral", "diff.same")
    ## plot.nodf(res.evo, "neutral", "diff.diff")
    ## mtext('Neutral', side=3, line=2, cex=1.5)
    plot.mod(res.evo, "qual", "same.same")
    plot.mod(res.evo, "qual", "same.diff")
    plot.mod(res.evo, "qual", "diff.same")
    plot.mod(res.evo, "qual", "diff.diff")
    plot.mod(res.evo, "quan", "same.same")
    plot.mod(res.evo, "quan", "same.diff")
    plot.mod(res.evo, "quan", "diff.same")
    plot.mod(res.evo, "quan", "diff.diff")
    mtext(xlabs, side=1, line=4, cex=1.5, adj=-2)
    ## plot.mod(res.evo, "neutral", "same.same")
    ## plot.mod(res.evo, "neutral", "same.diff")
    ## plot.mod(res.evo, "neutral", "diff.same")
    ## plot.mod(res.evo, "neutral", "diff.diff")
  }
  pdf.f(f, file= file.path(path.dir, sprintf("%s.pdf",
             paste(metric1, metric2, metric3,
                   sep=""))), width=7, height=6)
  
}
