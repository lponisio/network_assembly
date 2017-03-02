quartz()
layout(matrix(1:8,byrow=TRUE,nrow=2))

par(mar=c(2,4,3,0)+0.1)
plot.sim.results.bymat(sim.data, "same.same", "nodf", "fund", show.xlab=FALSE, show.ylab=TRUE)
mtext("Same-Same", 3)

par(mar=c(2,2,3,2)+0.1)
plot.sim.results.bymat(sim.data, "diff.same", "nodf", "fund", show.xlab=FALSE, show.ylab=FALSE)
mtext("Diff-same", 3)

par(mar=c(2,2,3,2)+0.1)
plot.sim.results.bymat(sim.data, "same.diff", "nodf", "fund", show.xlab=FALSE, show.ylab=FALSE)
mtext("Same-Diff", 3)

par(mar=c(2,0,3,4)+0.1)
plot.sim.results.bymat(sim.data, "diff.diff", "nodf", "fund", show.xlab=FALSE, show.ylab=FALSE)
mtext("diff-Diff", 3)

par(mar=c(4,4,1,0)+0.1)
plot.sim.results.bymat(sim.data, "same.same", "nodf", "eco", show.xlab=TRUE, show.ylab=TRUE)


par(mar=c(4,2,1,2)+0.1)
plot.sim.results.bymat(sim.data, "diff.same", "nodf", "eco", show.xlab=TRUE, show.ylab=FALSE)


par(mar=c(4,2,1,2)+0.1)
plot.sim.results.bymat(sim.data, "same.diff", "nodf", "eco", show.xlab=TRUE, show.ylab=FALSE)


par(mar=c(4,0,1,4)+0.1)
plot.sim.results.bymat(sim.data, "diff.diff", "nodf", "eco", show.xlab=TRUE, show.ylab=FALSE)


