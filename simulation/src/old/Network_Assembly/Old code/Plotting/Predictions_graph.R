source("/Users/laurenponisio/Documents/R General Functions/axis_arrows.R")

### by tree
cols <- c(hsv(0.4, s=1, alpha=0.1), hsv(0.5, s=1, alpha=0.1), hsv(0.7, s=1, alpha=0.1), hsv(0.8, s=1, alpha=0.1))

b.cols <- c(hsv(0.4, s=1), hsv(0.5, s=1), hsv(0.7, s=1), hsv(0.8, s=1))
	
xs <- c(2, 5, 0.5,0.5)

ys <-c(0.95, 0.5, 0.5, 0.05)


plot( xs, ys, xlim=c(0,6), ylim=c(-0.1,1.1), pch=21, cex=16, yaxt="n", xaxt="n", col=b.cols, bg=cols, ylab="Modularity", xlab="Nestedness")

legend("topright", legend=c("Tight coevo", "Diffuse coevo", "Geo-structuring", "Neutral"), pt.bg=cols, col=b.cols, pch=21, bty="n")

axis.arrows()





