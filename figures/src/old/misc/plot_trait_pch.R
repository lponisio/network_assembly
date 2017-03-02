plot.trait <- function(tree, traits, dir){
	standard.traits <- traits-min(traits)
	standard.traits <- standard.traits/max(standard.traits)
	plot(tree,show.tip.label=FALSE, direction=dir)
	tiplabels(pch=16,col="black",cex=standard.traits*3)
}

