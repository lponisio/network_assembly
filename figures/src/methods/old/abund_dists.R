library(RColorBrewer)
library(plotrix)
library(fields)
setwd("~/Dropbox/network_assembly")
cases <- expand.grid(mean.abund=log(seq(exp(0),exp(1), length=2)),
                     sd.abund = log(seq(exp(0),exp(1), length=2)))

setseed(10)
pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}
abund.plot <- function(cases, path){
  f <- function(){
    ## cases <- unique(cases[,2])
    abund.sim <- vector('list', length=nrow(cases))
    for(i in 1:nrow(cases)){
      abund.sim[[i]] <- rpoilog(10^5, mu=cases[i,1], sig=cases[i,2])
    }
    ## abund.sim <- abund.sim[sample(1:length(abund.sim), 3)]
    cols <- c('dodgerblue', 'dodgerblue4', 'darkolivegreen4',
              'darkolivegreen')
    ## plot(density(abund.sim[[1]]), col=cols[1], ylim=c(0,1.6),
    ##      xlim= c(0,20), xaxt='n', yaxt='n', ylab='', xlab='',
    ##      main='', type='l', lwd=2)
    plot(sort(abund.sim[[1]], TRUE), log='y', col=cols[1],
         xaxt='n', yaxt='n', ylab='', xlab='',
         ylim=c(1, 250),
         main='', type='l', lwd=2)  
    ## plot(hist(abund.sim[[1]], plot=FALSE,
    ##          breaks=seq(min(unlist(abund.sim)) -1,
    ##         max(unlist(abund.sim)) + 1, 1)), col=cols[1],
    ##      xaxt='n', yaxt='n', ylab='', xlab='',
    ##      main='', freq=FALSE)
    for(i in 2:length(abund.sim)){
      ## points(density(abund.sim[[i]]), type='l', col=cols[i+1],
      ## lwd=2)
      ## plot(hist(abund.sim[[i]], plot=FALSE,
      ##           breaks=seq(min(unlist(abund.sim)) -1,
      ##             max(unlist(abund.sim)) +1, 1)),
      ##      col=cols[i],
      ##      xaxt='n', yaxt='n', ylab='', xlab='',
      ##      main='', add=TRUE, freq=FALSE)
      points(sort(abund.sim[[i]], TRUE), type='l',
             col=cols[i],
             lwd=2)
    }
    mtext('Abundance (log)', 2, line=1, cex=1.5)
    mtext('Rank', 1, line=1, cex=1.5)

    ## gradient.rect(xleft=17, ybottom=0.95, xright=19, ytop= 1.5,
    ##               col=cols, gradient='y')
    
    ## text(labels=substitute(sigma == x, list(x=cases[length(cases)])),
    ##      x=15, y=1.50, cex=1.2)
    
    ## text(labels=substitute(sigma == x, list(x=cases[1])),
    ##      x=15, y=0.96, cex=1.2)

    rect(xleft=c(4*10^4, 4*10^4, 5*10^4, 5*10^4),
         xright=c(5*10^4, 5*10^4, 6*10^4, 6*10^4),
         ybottom=c(50, 100, 50, 100),
         ytop=c(100, 200, 100, 200), col=cols)
    
    text(labels=substitute(sigma == x, list(x=0)),
         x=4.5*10^4, y=250, cex=1.2)
    text(labels=substitute(sigma == x, list(x=1)),
         x=5.5*10^4, y=250, cex=1.2)
    text(labels=substitute(mu == x, list(x=1)),
         x=3.5*10^4, y=140, cex=1.2)
    text(labels=substitute(mu == x, list(x=0)),
         x=3.5*10^4, y=60, cex=1.2)
  }

  pdf.f(f, file= file.path(path, sprintf('%s.pdf', paste('abund',
             sep=''))), width=5, height=5)
  
}



abund.plot(cases, 'figures/methods')
