library(MASS)
library(ellipse)

CI95.2D <- function(x, y, col, b.col="black", n=50) {
  naOK <- !is.na(x) & !is.na(y)
  x <- x[naOK]
  y <- y[naOK]
  
  if(var(x) < 1e-04) {
    if(var(y) < 1e-04) {
      points(0, 0, bg=col, col= b.col, pch= 21)
    } else {
      yquant <- quantile(y, c(0.025, 0.975))
      ## rect(xleft=mean(x)-0.01*diff(par("usr")[1:2]),
      ##xright=mean(x)+0.01*diff(par("usr")[1:2]),
      ## ybottom=yquant[1],ytop=yquant[2],
      ## col=col,border=NA)

      elps <- ellipse(matrix(c(0.01*diff(par("usr")[1:2]), 0, 0,
                               quantile(y, 0.975) -
                               mean(y))/qnorm(0.975), 2,2)^2)
      
      elps[,1] <- elps[,1] + mean(x)
      elps[,2] <- elps[,2] + mean(y)
      
      polygon(elps, col= col, border= b.col)
    }
  } else if(var(y) < 1e-04) {
    xquant <- quantile(x,c(0.025, 0.975))
    ## rect(xleft=xquant[1],xright=xquant[2],
    ## ybottom=mean(y)-0.01*diff(par("usr")[3:4]),
    ## ytop=mean(y)+0.01*diff(par("usr")[3:4]),
    ## col=col,border=NA)
    
    elps <- ellipse(matrix(c(quantile(x,0.975) - mean(x),0,
                             0,0.01*diff(par("usr")[3:4]))/qnorm(0.975),
                           2,2)^2)
    
    elps[,1] <- elps[,1] + mean(x)
    elps[,2] <- elps[,2] + mean(y)
    
    polygon(elps,col=col,border=b.col)
  } else {
    xy.den <- kde2d(x, y, n=n)
    k <- sum(xy.den$z)
    
    zcut <- sapply(xy.den$z[xy.den$z < quantile(xy.den$z,0.75,
                                                na.rm=TRUE)],
                   function(x) sum(xy.den$z[xy.den$z > x])/k)
    
    this.z <- xy.den$z[xy.den$z < quantile(xy.den$z,0.75,
                                           na.rm=TRUE)][which.min(abs(zcut
                                             - 0.95))]
    
    browser()
    filled.contour(as.double(xy.den$x), as.double(xy.den$y), xy.den$z, 
                   zlim= as.double(c(this.z, 1.1*max(xy.den$z))),
                   #col = c(col, "transparent")
                   )
    
    contour(xy.den$x, xy.den$y, xy.den$z, levels= this.z, col=b.col,
            drawlabels=FALSE, add=TRUE)
  }
}
