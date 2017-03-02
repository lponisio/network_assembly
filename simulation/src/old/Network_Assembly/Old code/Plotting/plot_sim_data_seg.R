eg.phycor <- rep((7:4)/8,5)
eg.nodf <- eg.phycor + rep(1:5,rep(4,5))
eg.mat <- rep(c("fund","eco","samp1000","samp100"),5)


plot.sim.results.seg <- function(phy.cor,metric,mat,col,pchs=rep(1,4),...) {
	data2plot <- split(data.frame(phy=phy.cor,met=metric),mat)
	print(data2plot)
	names(pchs) <- c("fund","eco","samp1000","samp100")
	
	plot(phy.cor,metric,col=col,pch=pchs[mat],...)
	
	for(i in 1:3) {
		segments(data2plot[[i]]$phy,data2plot[[i]]$met,
				 data2plot[[i+1]]$phy,data2plot[[i+1]]$met,
				 col=col)
	}
	
}

plot.sim.results.seg <- function(phy.cor,metric,mat,col,pchs=rep(1,2),...) {
	data2plot <- split(data.frame(phy=phy.cor,met=metric),mat)
	names(pchs) <- c("fund","eco")
	
	plot(phy.cor,metric,col=col,pch=pchs[mat],...)
	
	for(i in 1:1) {
		segments(data2plot[[i]]$phy,data2plot[[i]]$met,
				 data2plot[[i+1]]$phy,data2plot[[i+1]]$met,
				 col=col)
	}
	
}

treetopos <- c("same.same", "diff.same", "same.diff", "diff.diff")

rules <- c( "trait.comp", "neutral","trait.narrow", "trait.wide", "trait.bar.comp","trait.bar.narrow", "trait.bar.wide", "trait.bar")
	
sim.data.abrev <- sim.data[sim.data$topo =="diff.same" &sim.data$age==100 & sim.data$npol == 30 & sim.data$mat == "eco" | sim.data$mat=="fund",]

layout(matrix(1:8, nrow=2))
for (i in 1:length(rules)){
	plot.sim.results.seg(sim.data.abrev$mean.cor[sim.data.abrev$link.rule==rules[i]],sim.data.abrev $nodf[sim.data.abrev$link.rule==rules[i]],sim.data.abrev$mat[sim.data.abrev$link.rule==rules[i]],col=hsv(0.6,alpha=0.1))
}
