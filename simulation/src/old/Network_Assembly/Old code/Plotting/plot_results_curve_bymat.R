

plot.sim.results.curve.byrule <- function(data, metric, cors, show.xlab=TRUE, show.ylab=TRUE, ...){
		
	rules <- c( "trait.comp", "neutral","trait.narrow", "trait.wide", "trait.bar.comp","trait.bar.narrow", "trait.bar.wide", "trait.bar")
	
	treetopo <- c("same.same", "diff.same", "same.diff", "diff.diff")
		
	cols <- brewer.pal(11, "Spectral")[8:11]
		
	data <- data[is.finite(data[,metric]),]
	
	layout(matrix(1:8, nrow=2))
	
	for (i in 1:length(rules)){
		for(j in 1:length(treetopos)){
			if(j==1){
			plot(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo[j], metric] ~
				data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo[j], cors],
				xlab=ifelse(show.xlab,"Interaction phylogenetic signal",""), ylab=ifelse(show.ylab,metric,""),
				col=cols[j], 
				ylim = range(data[, metric], na.rm=TRUE),
				xlim = range(data[,cors], na.rm=TRUE),
				type = "o", 
				main = rules[i],
				...)
			}else{	
				points(data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo[j], metric] ~ 
				data[data[,"link.rule"] == rules[i] & data[,"topo"] == treetopo[j], cors], 
				col=cols[j],
				type="o",
				...)
			}
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


