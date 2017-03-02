xy.sim.plot <- function(xsim,ysim,
					    xlim,ylim,
					    alpha=0.05,...) {
	if(missing(xlim)) xlim <- range(xsim)
	if(missing(ylim)) ylim <- range(ysim)
	
	mar <- par("mar")
	
	layout(matrix(c(2,1,0,3),2,2),widths=c(2,1.2),heights=c(1.2,2))
	
	par(mar=c(mar[1:2],0,0))
	plot(1,type="n",xlim=xlim,ylim=ylim,asp=1,...)
	for(i in 1:3) {
		xrng <- as.numeric(quantile(xsim[,i],c(alpha/2,1-alpha/2)))
		yrng <- as.numeric(quantile(ysim[,i],c(alpha/2,1-alpha/2)))
		rect(xrng[1],yrng[1],xrng[2],yrng[2],
			 col=hsv(0,0,c(0,0,0.5),alpha=0.7)[i],density=c(NA,9,NA)[i])
	}
	abline(0,1,lty=2)
	
	par(mar=c(0,mar[2:3],0))
	for(i in 1:3) {
		xden <- density(xsim[,i])
		if(i == 1) {
			plot(xden,col=hsv(0,0,c(0,0,0.5))[i],lwd=3,main="",xlab="",ylab="",   ########
				 xlim=xlim,axes=FALSE)    # ylim temporary fix only
		} else {
			lines(xden,col=hsv(0,0,c(0,0,0.5))[i],lwd=3)    ########
		}
		use.these <- which(xden$x >= quantile(xsim[,i],alpha/2) & xden$x <= quantile(xsim[,i],1-alpha/2))
			
		polygon(xden$x[c(min(use.these),use.these,max(use.these))],c(0,xden$y[use.these],0),
				col=hsv(0,0,c(0,0,0.5),alpha=0.7)[i],density=c(NA,9,NA)[i],border=NA)    ##########
		
		segments(x0=xden$x[range(use.these)],y0=c(0,0),y1=xden$y[range(use.these)],
				 col=hsv(0,0,c(0,0,0.5))[i])    ##############
	}
	
	par(mar=c(mar[1],0,0,mar[4]))
	for(i in 1:3) {
		yden <- density(ysim[,i])
		if(i == 1) {
			plot(yden$y,yden$x,type="l",col=hsv(0,0,c(0,0,0.5))[i],lwd=3,main="",xlab="",ylab="",   #########
				 ylim=ylim,axes=FALSE)
		}
		else {
			lines(yden$y,yden$x,col=hsv(0,0,c(0,0,0.5))[i],lwd=3)  ###############
		}
		
		use.these <- which(yden$x >= quantile(ysim[,i],alpha/2) & yden$x <= quantile(ysim[,i],1-alpha/2))
		
		polygon(y=yden$x[c(min(use.these),use.these,max(use.these))],x=c(0,yden$y[use.these],0),
				col=hsv(0,0,c(0,0,0.5),alpha=0.7)[i],density=c(NA,9,NA)[i],border=NA)  #########
		
		segments(y0=yden$x[range(use.these)],x0=c(0,0),x1=yden$y[range(use.these)],
				 col=hsv(0,0,c(0,0,0.5))[i])    ########
	}
	
	par(mar=mar)
}


##	this function expects the x-axis simulation and y-axis simulations
##	to be each a matrix having three columns (each column being a
##	different simulation)...something like:
xSim <- ySim <- t(matrix(rnorm(3*100),nrow=3)+c(0,10,100))

##	then you can plot that...but you might want to modify the
##	assumption of 3 simulations

xy.sim.plot(xSim,ySim)
