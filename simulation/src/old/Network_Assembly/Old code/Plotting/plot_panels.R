


treetopos <- c("same.same", "diff.same", "same.diff", "diff.diff")

matrixes <- c("fund", "eco", "samp100", "samp1000")

xlabs <- c(rep(FALSE,12), rep(TRUE,4))
ylabs <- c(rep(TRUE, 4), rep(FALSE,12))

library(RColorBrewer)

layout(matrix(1:16,nrow=4))

means1 <- means[means$age==100 & means$npol==30,]
meds1 <- meds[meds$age==100 & meds$npol==30,]
quartz()
for(i in 1:4){
	for(j in 1:4){
	plot.sim.results.bymat(means1, treetopos[i], "nodf", matrixes[j], "mean.cor", show.xlab=TRUE, show.ylab=TRUE)
	}
}


means2 <- means1[means1$mat != "samp",]
meds2 <- meds1[meds1$mat != "samp",]

means2$mat <- as.character(means2$mat)
means2$mat <- factor(means2$mat, levels = c("fund", "eco", "samp1000", "samp100"), ordered=TRUE)
means2 <- means2[order(means2$mat),]


meds2$mat <- as.character(meds2$mat)
meds2$mat <- factor(meds2$mat, levels = c("fund", "eco", "samp1000", "samp100"), ordered=TRUE)
meds2 <- meds2[order(meds2$mat),]


layout(matrix(1:4, nrow=1))

plot.sim.results.curve(means2, "same.same","nodf", "mean.cor", lwd="2", main="Same tree, same traits")
plot.sim.results.bymat(means2, "same.same","nodf", "fund", "mean.cor", pch=16, show.legend=TRUE)
plot.sim.results.bymat(means2, "same.same","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.same","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.same","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "diff.same","nodf", "mean.cor", lwd="2", main="Different tree, same traits")
plot.sim.results.bymat(means2, "diff.same","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "same.diff","nodf", "mean.cor", lwd="2", main="Same tree, different traits")
plot.sim.results.bymat(means2, "same.diff","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "diff.diff","nodf", "mean.cor", lwd="2", main="Different tree, different traits")
plot.sim.results.bymat(means2, "diff.diff","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "samp100","mean.cor", pch=16)


layout(matrix(1:4, nrow=1))

plot.sim.results.curve(means2, "same.same","nodf", "mean.cor", lwd="2", main="Same tree, same traits")
plot.sim.results.bymat(means2, "same.same","nodf", "fund", "mean.cor", pch=16, show.legend=TRUE)
plot.sim.results.bymat(means2, "same.same","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.same","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.same","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "diff.same","nodf", "mean.cor", lwd="2", main="Different tree, same traits")
plot.sim.results.bymat(means2, "diff.same","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.same","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "same.diff","nodf", "mean.cor", lwd="2", main="Same tree, different traits")
plot.sim.results.bymat(means2, "same.diff","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "same.diff","nodf", "samp100","mean.cor", pch=16)

plot.sim.results.curve(means2, "diff.diff","nodf", "mean.cor", lwd="2", main="Different tree, different traits")
plot.sim.results.bymat(means2, "diff.diff","nodf", "fund", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "eco", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "samp1000", "mean.cor", pch=16)
plot.sim.results.bymat(means2, "diff.diff","nodf", "samp100","mean.cor", pch=16)

## medians 
layout(matrix(1:4, nrow=1))

plot.sim.results.curve(meds2, "same.same","mod", "mean.cor", lwd="2", main="Same tree, same traits")


plot.sim.results.curve(meds2, "diff.same","mod", "mean.cor", lwd="2", main="Different tree, same traits")


plot.sim.results.curve(meds2, "same.diff","mod", "mean.cor", lwd="2", main="Same tree, different traits")


plot.sim.results.curve(meds2, "diff.diff","mod", "mean.cor", lwd="2", main="Different tree, different traits")




## by link rule

plot.sim.results.curve.byrule(means2, "nodf", "mean.cor", lwd=2)


legend("topright", legend=c("same.same", "diff.same", "same.diff", "diff.diff"), col=brewer.pal(11, "Spectral")[8:11], lwd=2, bty="n")


plot.sim.results.curve.byrule(means, "mod", "mean.cor", lwd=2)


legend("topright", legend=c("same.same", "diff.same", "same.diff", "diff.diff"), col=brewer.pal(11, "Spectral")[8:11], lwd=2, bty="n")
		
		
plot.sim.results.curve.byrule(meds2, "nodf", "mean.cor", lwd=2)


legend("topright", legend=c("same.same", "diff.same", "same.diff", "diff.diff"), col=brewer.pal(11, "Spectral")[8:11], lwd=2, bty="n")


plot.sim.results.curve.byrule(meds, "mod", "mean.cor", lwd=2)


legend("topright", legend=c("same.same", "diff.same", "same.diff", "diff.diff"), col=brewer.pal(11, "Spectral")[8:11], lwd=2, bty="n")
		

		

