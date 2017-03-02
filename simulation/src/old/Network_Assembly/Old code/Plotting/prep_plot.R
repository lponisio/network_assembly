library(RColorBrewer)

prep.plot <- function(simres, metric){
	browser()
	means.sim <- aggregate(list(metrics=simres[[1]][,metric]), 
	by= list(topo=simres[[2]][,"topo"],
          link.rule=simres[[2]][,"link.rule"],
          mat=simres[[2]][,"mat"], mu=unlist(simres[[2]][,"mu"]),
          lambda=unlist(simres[[2]][,"lambda"])), 
		mean, na.rm=TRUE)

	cats <- paste((means.sim$topo), (means.sim$link.rule), (means.sim$mat))

	split.means <- split(means.sim, cats)

	prep.mats <- lapply(split.means, FUN=function(X){
		test <- matrix(X[,"metrics"], nrow=10, byrow=TRUE)
		rownames(test) <-unique(X[,"lambda"])
		colnames(test) <- unique(X[,"mu"])
		return(test)
	})

	funds <- prep.mats[grep("fund", names(prep.mats))]
	ecos <- prep.mats[grep("eco", names(prep.mats))]

	par(mfrow=c(4,3),las=1,cex.axis=1,cex.lab=1)
	par(mar=c(2,3,2,2))
	
	n.breaks <- 20
	
	all.dats <- rbind(unlist(ecos), unlist(funds))
	all.dats <- all.dats[is.finite(all.dats)]

	lapply(funds, image, ylab="mu", xlab="lambda", col=rainbow(n.breaks-1, end=0.85), breaks=seq(from=0, to=max(all.dats),length=n.breaks))
	
	quartz()
	
	par(mfrow=c(4,3),las=1,cex.axis=1,cex.lab=1)
	par(mar=c(2,3,2,2))
	
	lapply(ecos, image, ylab="mu", xlab="lambda", col=rainbow(n.breaks-1, end=0.85), breaks=seq(from=0, to=max(all.dats),length=n.breaks))
	
	return(list(funds=funds, ecos=ecos))
}
