##plot by topology combinations
##xdatas should be a list with each element equal to the fund data for each tree combination, same with the ys

cor.plots <- function(xdatas, ydata1s, ydata2s, ydata3s, trait){ 
	cols <- c(hsv(0.5, s=0.5, alpha=0.99), hsv(0.1, s=0.5, alpha=0.99), hsv(1, s=0.5, alpha=0.99))
	
	layout(matrix(1:(length(trait)*length(xdatas)), ncol=length(trait), nrow=length(xdatas),byrow=TRUE))	
	
	mains <- c("Complimentary Trait", "Barrior Trait", "Neutral Trait")
	
	old.mar <- par("mar")
	
	
	for(i in 1:length(xdatas)){
	
		xdata1 <- xdatas[[i]]
		ydata1 <- ydata1s[[i]]
		ydata2 <- ydata2s[[i]]
		ydata3 <- ydata3s[[i]]
	
		allys <- rbind(ydata1, ydata2, ydata3)
		
	###first trait
		for(j in 1:length(trait)){
			
			if(i==1){
			par(mar = old.mar + c(-2.5,0,-2, -2))
			}
			else if(j == 1 & i != 4){
			par(mar = old.mar + c(-2,0,-2.5,-2))
			}
			else if(i==4){
			par(mar = old.mar + c(-1,0,-3.5,-2))
			}
			else{
			par(mar = old.mar + c(-2,0,-2.5,-2))
			}
		
			if(i==1){
				
			}
			
			plot(xdata1$nodf, ydata1$nodf, type="n", ylim=c(0,100), xlim=c(0,100), xlab="", ylab="", xaxt="n", yaxt="n"); abline(0,1)
			
			if(i==1){
				mtext(mains[[j]], 3, line=0.5)
			}
			
			if(i == length(xdatas)){
				axis(side=1)
				mtext("Nestedness fundamental", side=1, line=2.5, cex=0.55)
			}
			
			if(j == 1){
				axis(side=2)
				mtext("Nestedness ecological and sampling", side=2.5, line=2, cex=0.55)
			}
	
			CI95.2D(xdata1$nodf[xdata1$link.rule == trait[j]],
				ydata2$nodf[ydata2$link.rule == trait[j]],
				col=cols[2])
		
			CI95.2D(xdata1$nodf[xdata1$link.rule == trait[j]],
				ydata3$nodf[ydata3$link.rule == trait[j]],
				col=cols[3])
		
			CI95.2D(xdata1$nodf[xdata1$link.rule == trait[j]],
				ydata1$nodf[ydata1$link.rule == trait[j]],
				col=cols[1])
	
			if(i == 1 & j == 3){
				legend("topleft", legend=c("Ecological", "Good sampling", "Poor sampling"), col=cols, pch=16, bty="n")	
			}	##close if
		}	## close j loop 
	} ## close i loop	
	
	quartz()
	layout(matrix(1:(length(trait)*length(xdatas)), ncol=length(trait), nrow=length(xdatas),byrow=TRUE))
	
	for(i in 1:length(xdatas)){
	
		xdata1 <- xdatas[[i]]
		ydata1 <- ydata1s[[i]]
		ydata2 <- ydata2s[[i]]
		ydata3 <- ydata3s[[i]]
		
		xdata1$mod[xdata1$mod < 0] <- 0
		ydata1$mod[ydata1$mod < 0] <- 0
		ydata2$mod[ydata2$mod < 0] <- 0
		ydata3$mod[ydata3$mod < 0] <- 0
	
		allys <- rbind(ydata1, ydata2, ydata3)		
	
		for(j in 1:length(trait)){
		if(i==1){
			par(mar = old.mar + c(-2.5,0,-2, -2))
			}
			else if(j == 1 & i != 4){
			par(mar = old.mar + c(-2,0,-2.5,-2))
			}
			else if(i==4){
			par(mar = old.mar + c(-1,0,-3.5,-2))
			}
			else{
			par(mar = old.mar + c(-2,0,-2.5,-2))
			}
		
			if(i==1){
				
			}
			
			plot(xdata1$mod, ydata1$mod, type="n", ylim=c(0,1), xlim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n"); abline(0,1)
			
			if(i==1){
				mtext(mains[[j]], 3, line=0.5)
			}
			
			if(i == length(xdatas)){
				axis(side=1)
				mtext("Modularity fundamental", side=1, line=2.5, cex=0.55)
			}
			
			if(j == 1){
				axis(side=2)
				mtext("Modularity ecological and sampling", side=2.5, line=2, cex=0.55)
			}
	
		CI95.2D(xdata1$mod[xdata1$link.rule == trait[j]],
			ydata2$mod[ydata2$link.rule == trait[j]],
			col=cols[2])
		
		CI95.2D(xdata1$mod[xdata1$link.rule == trait[j]],
			ydata3$mod[ydata3$link.rule == trait[j]],
			col=cols[3])
		
		CI95.2D(xdata1$mod[xdata1$link.rule == trait[j]],
			ydata1$mod[ydata1$link.rule == trait[j]],
			col=cols[1])
		
		if(i == 1 & j == 3){
			legend("topleft", legend=c("Ecological", "Good sampling", "Poor sampling"), col=cols, pch=16, bty="n")	
		}	
				
		
		}
	}	
}


