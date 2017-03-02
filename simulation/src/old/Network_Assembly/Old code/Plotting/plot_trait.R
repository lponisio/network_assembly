plot.trait <- function(tree, traits, dir){
	cInt <- classIntervals(traits,24,style="fisher")
	cPal <- tim.colors(24)
	#plot(cInt,cPal)
	these.col <- findColours(cInt,cPal)

	plot(tree,show.tip.label=FALSE, direction=dir)
	tiplabels(pch=16,col=these.col,cex=1.5)

}

