
library(RColorBrewer)


plot.sim.results.curve <- function(data, treetopo, metric, cors, show.xlab=TRUE, show.ylab=TRUE, ...){
		
	rules <- c( "trait.comp", "neutral","trait.narrow", "trait.wide", "trait.bar.comp","trait.bar.narrow", "trait.bar.wide", "trait.bar")
		
	cols <- brewer.pal(8, "Dark2")
		
	data <- data[is.finite(data[,metric]),]
	
	for (i in 1:length(rules)){
		if(i ==1){
		plot(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo, metric] ~
			data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo, cors],
			xlab=ifelse(show.xlab,"Interaction phylogenetic signal",""), ylab=ifelse(show.ylab,metric,""),
			col=cols[i], 
			ylim = range(data[, metric], na.rm=TRUE),
			xlim = range(data[,cors], na.rm=TRUE),
			type= "l",
			...)
		} else{
			points(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo, metric] ~ 
			data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo, cors], 
			col=cols[i],
			type="l",
			...)
		}
	
	}
}




library(RColorBrewer)


plot.sim.results.bymat <- function(data, treetopo, metric, matrix.type, cors, show.xlab=TRUE, show.ylab=TRUE, show.legend=FALSE, ...){
		
	rules <- c( "trait.comp", "neutral","trait.narrow", "trait.wide", "trait.bar.comp","trait.bar.narrow", "trait.bar.wide", "trait.bar")
		
	cols <- brewer.pal(8, "Dark2")
		
	data <- data[is.finite(data[,metric]),]
	
	for (i in 1:length(rules)){
	
		points(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, cors], 
		col=cols[i],...)
	}
	if(show.legend==TRUE){
		legend("topright", legend= rules, col=cols, lwd=2, bty="n")
	}
}


