##only fund vs. other matrices plaotted, not eco vs others


cor.plots <- function(xdata1, ydata1, ydata2, ydata3, trait){
	
	cols <- c(hsv(0.5, alph=0.75), hsv(0.6, alph=0.75), hsv(0.8, alph=0.75))
	allxs <- xdata1
	allys <- rbind(ydata1, ydata2, ydata3)
	
	layout(matrix(1:3, ncol=3, nrow=1))
	
	
	###first trait
	plot(xdata1$nodf, ydata1$nodf, type="n", ylim=c(0,100), xlim=c(0,100), xlab="Nestedness Fundamental", ylab="Nestedness ecological and sampling", main="Complimentary Trait")
	abline(0,1)
	
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[1]],
		ydata2$nodf[ydata2$link.rule == trait[1]],
		col=cols[2])
		
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[1]],
		ydata3$nodf[ydata3$link.rule == trait[1]],
		col=cols[3])
		
		CI95.2D(xdata1$nodf[xdata1$link.rule == trait[1]],
		ydata1$nodf[ydata1$link.rule == trait[1]],
		col=cols[1])
	legend("topleft", legend=c("Ecological", "Good sampling", "Poor sampling"), col=cols, pch=16, bty="n")	
		
			
		##second trait
		plot(xdata1$nodf, ydata1$nodf, type="n", ylim=c(0,100), xlim=c(0,100), xlab="Nestedness Fundamental", ylab="Nestedness ecological and sampling", main="Barrior Trait")
	abline(0,1)
	
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[2]],
		ydata2$nodf[ydata2$link.rule == trait[2]],
		col=cols[2])
		
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[2]],
		ydata3$nodf[ydata3$link.rule == trait[2]],
		col=cols[3])
		
		CI95.2D(xdata1$nodf[xdata1$link.rule == trait[2]],
		ydata1$nodf[ydata1$link.rule == trait[2]],
		col=cols[1])
	
		
		##third trait
		
				
		plot(xdata1$nodf, ydata1$nodf, type="n", ylim=c(0,100), xlim=c(0,100), xlab="Nestedness Fundamental", ylab="Nestedness ecological and sampling",main="Neutral Trait")
	abline(0,1)
	
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[3]],
		ydata2$nodf[ydata2$link.rule == trait[3]],
		col=cols[2])
		
	CI95.2D(xdata1$nodf[xdata1$link.rule == trait[3]],
		ydata3$nodf[ydata3$link.rule == trait[3]],
		col=cols[3])
		
		CI95.2D(xdata1$nodf[xdata1$link.rule == trait[3]],
		ydata1$nodf[ydata1$link.rule == trait[3]],
		col=cols[1])
	
		
	quartz()
			
	layout(matrix(1:3, ncol=3, nrow=1))
	
	plot(xdata1$mod, ydata1$mod, type="n", ylim=c(0,1), xlim=c(0,1), xlab="Modularity Fundamental", ylab="Modularity ecological and sampling", main= "Complimentary Trait")
	abline(0,1)
	
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[1]],
		ydata2$mod[ydata2$link.rule == trait[1]],
		col=cols[2])
		
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[1]],
		ydata3$mod[ydata3$link.rule == trait[1]],
		col=cols[3])
		
		CI95.2D(xdata1$mod[xdata1$link.rule == trait[1]],
		ydata1$mod[ydata1$link.rule == trait[1]],
		col=cols[1])
	legend("topleft", legend=c("Ecological", "Good sampling", "Poor sampling"), col=cols, pch=16, bty="n")	
				
		
	##second trait
	
	xdata1$mod[xdata1$mod < 0] <- 0
	ydata1$mod[ydata1$mod < 0] <- 0
	ydata2$mod[ydata2$mod < 0] <- 0
	ydata3$mod[ydata3$mod < 0] <- 0
	
	
	plot(xdata1$mod, ydata1$mod, type="n", ylim=c(0,1), xlim=c(0,1), xlab="Modularity Fundamental", ylab="Modularity ecological and sampling", main="Barrior Trait")
	abline(0,1)
	
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[2]],
		ydata2$mod[ydata2$link.rule == trait[2]],
		col=cols[2])
		
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[2]],
		ydata3$mod[ydata3$link.rule == trait[2]],
		col=cols[3])
		
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[2]],
		ydata1$mod[ydata1$link.rule == trait[2]],
		col=cols[1])
	
		##third trait
		
	plot(xdata1$mod, ydata1$mod, type="n", ylim=c(0,1), xlim=c(0,1), xlab="Modularity Fundamental", ylab="Modularity ecological sampling", main="Neutral Trait")
	abline(0,1)
				
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[3]],
		ydata2$mod[ydata2$link.rule == trait[3]],
		col=cols[2])
		
	CI95.2D(xdata1$mod[xdata1$link.rule == trait[3]],
		ydata3$mod[ydata3$link.rule == trait[3]],
		col=cols[3])
		
		CI95.2D(xdata1$mod[xdata1$link.rule == trait[3]],
		ydata1$mod[ydata1$link.rule == trait[3]],
		col=cols[1])
		
		
}


cor.plots(fund.dd, eco.dd, samp1000.dd, samp100.dd, c("trait.narrow", "trait.bar", "neutral"))
cor.plots(fund.ss, eco.ss, samp1000.ss, samp100.ss, c("trait.narrow", "trait.bar", "neutral"))
cor.plots(fund.ds, eco.ds, samp1000.ds, samp100.ds, c("trait.narrow", "trait.bar", "neutral"))
cor.plots(fund.sd, eco.sd, samp1000.sd, samp100.sd, c("trait.narrow", "trait.bar", "neutral"))

