##plot modulaity in y and nodf on x

modVnodf.plots.bymat <- function(dat.ss, dat.ds, dat.sd, dat.dd, trait){ ##ssamental, dslogical, good sampling, poor sampling
	
	cols <- c(hsv(0.4, s=1, alpha=0.1), hsv(0.5, s=1, alpha=0.1), hsv(0.7, s=1, alpha=0.1), hsv(0.8, s=1, alpha=0.1))
	b.cols <- c(hsv(0.4, s=1), hsv(0.5, s=1), hsv(0.7, s=1), hsv(0.8, s=1))
	
	layout(matrix(1:(length(trait)*length(dat.ss)), ncol=length(trait), nrow=length(dat.ss),byrow=TRUE))	
	
	mains <- c("Complimentary Trait", "Barrier Trait", "Neutral Trait")
	
	old.mar <- par("mar")
	
	
	for(i in 1:length(dat.ss)){
	
		ss <- dat.ss[[i]]
		ds <- dat.ds[[i]]
		sd <- dat.sd[[i]]
		dd <- dat.dd[[i]]
	
		ss$mod[ss$mod < 0] <- 0
		ds$mod[ds$mod < 0] <- 0
		sd$mod[sd$mod < 0] <- 0
		dd$mod[dd$mod < 0] <- 0
	
		
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
			
			plot(ss$nodf, ss$mod, type="n", ylim=c(0,1), xlim=c(0,100), xlab="", ylab="", xaxt="n", yaxt="n")
			
			if(i==1){
				mtext(mains[[j]], 3, line=0.5)
			}
			
			if(i == length(dat.ss)){
				axis(side=1)
				mtext("Nestedness", side=1, line=2.5, cex=0.55)
			}
			
			if(j == 1){
				axis(side=2)
				mtext("Modularity", side=2, line=2.5, cex=0.55)
			}
	
			CI95.2D(ss$nodf[ss$link.rule == trait[j]],
				ss$mod[ss$link.rule == trait[j]],
				col=cols[1], b.col=b.cols[1])
		
			CI95.2D(ds$nodf[ds$link.rule == trait[j]],
				ds$mod[ds$link.rule == trait[j]],
				col=cols[2], b.col=b.cols[2])
		
			CI95.2D(sd$nodf[sd$link.rule == trait[j]],
				sd$mod[sd$link.rule == trait[j]],
				col=cols[3], b.col=b.cols[3])
				
			CI95.2D(dd$nodf[dd$link.rule == trait[j]],
				dd$mod[dd$link.rule == trait[j]],
				col=cols[4], b.col=b.cols[4])	
	
			if(i == 1 & j == 3){
				legend("topleft", legend=c("Tight coeco", "Diffuse coeco", "Geo-structuring", "Neutral"), pt.bg=cols, col=b.cols, pch=21, bty="n")	
			}	##close if
		}	## close j loop 
	} ## close i loop	
	
}


