


library(RColorBrewer)


plot.sim.results.bymat <- function(data, treetopo, metric, matrix.type, show.xlab=TRUE, show.ylab=TRUE, ...){
		
	cols <- brewer.pal(8, "Dark2")
		
	data <- data[is.finite(data[,metric]),]
	
	##narrow
	points(data[data[,"link.rule"] == "trait.narrow" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~
		data[data[,"link.rule"] == "trait.narrow" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"],
		...)
	
	#comp 		
	points(data[data[,"link.rule"] == "trait.comp" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.comp" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[2],...)
		
	#neutral	
	if(matrix.type != "fund") {
		points(data[data[,"link.rule"] == "neutral" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
			data[data[,"link.rule"] == "neutral" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
			col=cols[3],...)
	}
	
	#bar
	points(data[data[,"link.rule"] == "trait.bar" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.bar" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[4],...)
	
	#wide
	points(data[data[,"link.rule"] == "trait.wide" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.wide" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[5],...)
	
	##bar/wide	
	points(data[data[,"link.rule"] == "trait.bar.wide" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.bar.wide" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[6],...)
		
	##bar.narrow	
	points(data[data[,"link.rule"] == "trait.bar.narrow" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.bar.narrow" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[7],...)
	
	##bar/comp	
	points(data[data[,"link.rule"] == "trait.bar.comp" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, metric] ~ 
		data[data[,"link.rule"] == "trait.bar.comp" & data[,"topo"] == treetopo & data[,"mat"] == matrix.type, "mean.cor"], 
		col=cols[8],...)
	
}

