
library(RColorBrewer)


plot.sim.results.bymat <- function(data, treetopo, metric, matrix.type, cors, show.xlab=TRUE, show.ylab=TRUE, ...){
		
	rules <- c( "trait.comp", "neutral","trait.narrow", "trait.wide", "trait.bar.comp","trait.bar.narrow", "trait.bar.wide", "trait.bar")
		
	cols <- brewer.pal(8, "Dark2")
		
	data <- data[is.finite(data[,metric]),]
	
	for (i in 1:length(rules)){
		if(i ==1){
		plot(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~
			data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, cors],
			xlab=ifelse(show.xlab,"Pollinator phylogenetic signal",""), ylab=ifelse(show.ylab,metric,""),
			col=cols[i], 
			ylim = range(data[, metric], na.rm=TRUE),
			xlim = range(data[,cors], na.rm=TRUE),
			...)
		} else{
			points(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
			data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, cors], 
			col=cols[i],...)
		}
	
	}
}
