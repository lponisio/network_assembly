library(RColorBrewer)

plot.rho <- function(cors.nodf, cors.mod){
	
	cols <- brewer.pal(11, "Spectral")[8:11]
	
	layout(matrix(1:6, nrow=3, ncol=2))
	old.mar <- par("mar")
	par(mar = old.mar + c(-2.5,0,-2, -2))
	
	
	for (i in 1:3){
		
		for (j in 1:4){
			
			cors.nodf[[j]][!is.finite(cors.nodf[[j]])] <- 0
			
			if(j == 1){
				plot(cors.nodf[[j]][i, c(1,4,5)], xlim=c(1,3), type="o", ylab="cor", ylim=c(-1,1), col=cols[j], lwd=2, cex=2, xaxt="n", xlab="")
			} else {
				points(cors.nodf[[j]][i, c(1,4,5)], xlim=c(1,3), type="o", ylab="cor", col=cols[j],lwd=2, cex=2)
			}
			
			if(i == 1){
				mtext("Nestedness")
			}
			if(i == 3){
				axis(side=1, labels=c("Fund", "Good samp", "Poor samp"), at=c(1,2,3))
			}
			
		} #close j
	} #close i
	
	par(mar = old.mar + c(-2.5,-2,-2, 0))
	for (i in 1:3){
		
		for (j in 1:4){
			
			cors.mod[[j]][!is.finite(cors.mod[[j]])] <- 0
			
			if(j == 1){
				plot(cors.mod[[j]][i, c(1,4,5)], xlim=c(1,3), type="o", ylim=c(-1,1), col=cols[j], lwd=2, cex=2, xaxt="n", xlab="", ylab="", yaxt="n")
			} else {
				points(cors.mod[[j]][i, c(1,4,5)], xlim=c(1,3), type="o", ylab="cor", col=cols[j],lwd=2, cex=2)
			}
			
			if(i == 1){
				mtext("Modularity")
			}
			if(i == 3){
				axis(side=1, labels=c("Fund", "Good samp", "Poor samp"), at=c(1,2,3))
			}
			
		} #close j
	} #close i
	
	legend("topright", legend=c("Tight coeco", "Diffuse coeco", "Geo-structuring", "Neutral"), col=cols, lwd=2, bty="n")
}

