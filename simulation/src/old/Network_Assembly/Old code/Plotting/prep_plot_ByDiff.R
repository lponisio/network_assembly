##takes the simulation results in list form where the first element of the list if the data and the second is the catagories
## calculates mu - lambda and plots the mean metric by the diff

prep.plot.by.diff <- function(simres, metric, legend.loc, place.legend){
	
	simres[[2]] <- cbind(simres[[2]], unlist(simres[[2]][,"mu"]) - unlist(simres[[2]][,"lambda"]))
	names(simres[[2]]) <- c(names(simres[[2]][-7]), "diff")

	means.sim <- aggregate(list(metrics=simres[[1]][,metric]), 
	by=list(topo=simres[[2]][,"topo"],link.rule=simres[[2]][,"link.rule"],mat=simres[[2]][,"mat"], diff=simres[[2]][,"diff"]), 
		mean, na.rm=TRUE)
	
	split.mats <- split(means.sim,  means.sim$mat)
	
	split.means <- lapply(split.mats, FUN=function(x){
		cats <- paste((x$topo), (x$link.rule))
		split(x, cats)
	})
	
	cols <-brewer.pal(4, "RdYlBu")
	quartz(8,10)
	layout(matrix(1:2,nrow=1))
	
	for(j in 1:length(split.means)){
		for(i in 1:length(split.means[[j]])){
			
			# print(split.means[[j]][[i]])
			col <- switch(as.character(split.means[[j]][[i]]$topo[1]),
									same.same=cols[1],
									same.diff=cols[2],
									diff.diff=cols[3],
									diff.same=cols[4])
			lty <- switch(as.character(split.means[[j]][[i]]$link.rule[1]),
									trait.bar=2,
									trait.narrow=1,
									neutral=2)
			if(i == 1){
				plot(y=split.means[[j]][[i]][,"metrics"], x=split.means[[j]][[i]][,"diff"], ylim=c(0,max(means.sim$metrics[is.finite(means.sim$metrics)])), xlim=c(min(simres[[2]][,"diff"]),max(simres[[2]][,"diff"])), type="l",
						col=col,lty=lty,main=split.means[[j]][[i]]$mat[1], lwd=2, ylab=metric, xlab="Mu - lamdba")
				if(j==1 & i==1 & place.legend == "first"){
						legend(legend.loc, legend=c("Coevolution and cospeciation", "Cospeciation only", "Neutral", "Coevolution only", "Barrier trait", "Complimentary trait", "Neutral trait"), col=c(cols, "black", "black", "black"), lty=c(rep(1,4), 2,1,3), bty="n" )
				}		
			} else{
				points(y=split.means[[j]][[i]][,"metrics"], x=split.means[[j]][[i]][,"diff"], type="l",
							col=col,lty=lty, lwd=2)
			} #close else
		} #close j
	} #close i
	if( place.legend == "second"){
		legend(legend.loc, legend=c("Coevolution and cospeciation", "Cospeciation only", "Neutral", "Coevolution only", "Barrier trait", "Complimentary trait", "Neutral trait"), col=c(cols, "black", "black", "black"), lty=c(rep(1,4), 2,1,3), bty="n" )
	} # close if
	return(split.means)
} #close function