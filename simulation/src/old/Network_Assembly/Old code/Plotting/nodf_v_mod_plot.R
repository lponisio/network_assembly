##plot modulaity in y and nodf on x

modVnodf.plots <- function(funds, ecos, samp1000s, samp100s, trait){ ##fundamental, ecological, good sampling, poor sampling
	
	cols <- c(hsv(0.3, s=1, alpha=0.1), hsv(0.5, s=1, alpha=0.1), hsv(0.1, s=1, alpha=0.1), hsv(1, s=1, alpha=0.1))
	b.cols <- c(hsv(0.3, s=1), hsv(0.5, s=1), hsv(0.1, s=1), hsv(1, s=1))
	layout(matrix(1:(length(trait)*length(funds)), ncol=length(trait), nrow=length(funds),byrow=TRUE))	
	
	mains <- c("Complimentary Trait", "Barrier Trait", "Neutral Trait")
	
	old.mar <- par("mar")
	
	
	for(i in 1:length(funds)){
	
		fund <- funds[[i]]
		eco <- ecos[[i]]
		samp1000 <- samp1000s[[i]]
		samp100 <- samp100s[[i]]
	
		fund$mod[fund$mod < 0] <- 0
		eco$mod[eco$mod < 0] <- 0
		samp1000$mod[samp1000$mod < 0] <- 0
		samp100$mod[samp100$mod < 0] <- 0
	
		
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
			
			plot(eco$nodf, eco$mod, type="n", ylim=c(0,1), xlim=c(0,100), xlab="", ylab="", xaxt="n", yaxt="n")
			
			if(i==1){
				mtext(mains[[j]], 3, line=0.5)
			}
			
			if(i == length(funds)){
				axis(side=1)
				mtext("Nestedness", side=1, line=2.5, cex=0.55)
			}
			
			if(j == 1){
				axis(side=2)
				mtext("Modularity", side=2, line=2.5, cex=0.55)
			}
	
			CI95.2D(fund$nodf[fund$link.rule == trait[j]],
				fund$mod[fund$link.rule == trait[j]],
				col=cols[1], b.col=b.cols[1])
		
			CI95.2D(eco$nodf[eco$link.rule == trait[j]],
				eco$mod[eco$link.rule == trait[j]],
				col=cols[2], b.col=b.cols[2])
		
			CI95.2D(samp1000$nodf[samp1000$link.rule == trait[j]],
				samp1000$mod[samp1000$link.rule == trait[j]],
				col=cols[3], b.col=b.cols[3])
				
			CI95.2D(samp100$nodf[samp100$link.rule == trait[j]],
				samp100$mod[samp100$link.rule == trait[j]],
				col=cols[4], b.col=b.cols[4])	
	
			if(i == 1 & j == 3){
				legend("topleft", legend=c("Fundamental", "Ecological", "Good sampling", "Poor sampling"), pt.bg=cols, col=b.cols, pch=21, bty="n")	
			}	##close if
		}	## close j loop 
	} ## close i loop	
	
}


