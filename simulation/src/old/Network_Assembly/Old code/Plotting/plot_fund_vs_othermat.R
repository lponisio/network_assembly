load("loaded_data_prelim_4.RData")
source("/Users/laurenponisio/Documents/Network_Assembly/CI95_2D.R")
abrev.data <- sim.data.data[sim.data.cat$age == 100 & sim.data.cat$npol == 30,]

abrev.cat <- sim.data.cat[sim.data.cat$age == 100 & sim.data.cat$npol == 30,]

abrev.all <- as.data.frame(abrev.data)

abrev.all <- cbind(abrev.all, abrev.cat)

## fundamentals 
fund.ss <- abrev.all[abrev.all$mat == "fund" & abrev.all$topo == "same.same",]

fund.sd <- abrev.all[abrev.all$mat == "fund" & abrev.all$topo == "same.diff",]

fund.dd <- abrev.all[abrev.all$mat == "fund" & abrev.all$topo == "diff.diff",]

fund.ds <- abrev.all[abrev.all$mat == "fund" & abrev.all$topo == "diff.same",]

## ecological 
eco.ss <- abrev.all[abrev.all$mat == "eco" & abrev.all$topo == "same.same",]

eco.sd <- abrev.all[abrev.all$mat == "eco" & abrev.all$topo == "same.diff",]

eco.dd <- abrev.all[abrev.all$mat == "eco" & abrev.all$topo == "diff.diff",]

eco.ds <- abrev.all[abrev.all$mat == "eco" & abrev.all$topo == "diff.same",]

##sampling 100

samp100.ss <- abrev.all[abrev.all$mat == "samp100" & abrev.all$topo == "same.same",]

samp100.sd <- abrev.all[abrev.all$mat == "samp100" & abrev.all$topo == "same.diff",]

samp100.dd <- abrev.all[abrev.all$mat == "samp100" & abrev.all$topo == "diff.diff",]

samp100.ds <- abrev.all[abrev.all$mat == "samp100" & abrev.all$topo == "diff.same",]

##sampling 1000

samp1000.ss <- abrev.all[abrev.all$mat == "samp1000" & abrev.all$topo == "same.same",]

samp1000.sd <- abrev.all[abrev.all$mat == "samp1000" & abrev.all$topo == "same.diff",]

samp1000.dd <- abrev.all[abrev.all$mat == "samp1000" & abrev.all$topo == "diff.diff",]

samp1000.ds <- abrev.all[abrev.all$mat == "samp1000" & abrev.all$topo == "diff.same",]

###

plot(fund.dd[,"nodf"], samp1000.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "neutral"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "trait.narrow"],
		col=hsv(0.5, alpha=0.5))
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.bar"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "trait.bar"],
		col=hsv(1, alpha=0.5))
	
						

plot(fund.dd[,"nodf"], samp100.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp100.dd$nodf[samp100.dd$link.rule == "neutral"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.dd$nodf[samp100.dd$link.rule == "trait.narrow"],
		col=hsv(0.5, alpha=0.5))
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.bar"],
		samp100.dd$nodf[samp100.dd$link.rule == "trait.bar"],
		col=hsv(1, alpha=0.5))
	


plot(fund.dd[,"nodf"], eco.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		eco.dd$nodf[eco.dd$link.rule == "neutral"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.dd$nodf[eco.dd$link.rule == "trait.narrow"],
		col=hsv(0.5, alpha=0.5))
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.bar"],
		eco.dd$nodf[eco.dd$link.rule == "trait.bar"],
		col=hsv(1, alpha=0.5))
	





plot(fund.dd[,"nodf"], samp1000.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "neutral"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		eco.dd$nodf[eco.dd$link.rule == "neutral"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp100.dd$nodf[samp100.dd$link.rule == "neutral"],
		col=hsv(0.9, alpha=0.5))



plot(fund.dd[,"nodf"], samp1000.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "neutral"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		eco.dd$nodf[eco.dd$link.rule == "neutral"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "neutral"],
		samp100.dd$nodf[samp100.dd$link.rule == "neutral"],
		col=hsv(0.9, alpha=0.5))
		
plot(fund.dd[,"nodf"], samp1000.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "trait.narrow"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.dd$nodf[eco.dd$link.rule == "trait.narrow"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.dd$nodf[samp100.dd$link.rule == "trait.narrow"],
		col=hsv(0.9, alpha=0.5))
		
plot(fund.dd[,"nodf"], samp1000.dd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp1000.dd$nodf[samp1000.dd$link.rule == "trait.narrow"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.dd$nodf[eco.dd$link.rule == "trait.narrow"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.dd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.dd$nodf[samp100.dd$link.rule == "trait.narrow"],
		col=hsv(0.9, alpha=0.5))
		
plot(fund.ss[,"nodf"], samp1000.ss[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.ss$nodf[fund.ss$link.rule == "trait.narrow"],
		samp1000.ss$nodf[samp1000.ss$link.rule == "trait.narrow"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.ss$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.ss$nodf[eco.ss$link.rule == "trait.narrow"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.ss$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.ss$nodf[samp100.ss$link.rule == "trait.narrow"],
		col=hsv(0.9, alpha=0.5))
		


plot(fund.sd[,"nodf"], samp1000.sd[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.sd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp1000.sd$nodf[samp1000.dd$link.rule == "trait.narrow"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.sd$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.sd$nodf[eco.dd$link.rule == "trait.narrow"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.sd$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.sd$nodf[samp100.dd$link.rule == "trait.narrow"],
		col=hsv(0.9, alpha=0.5))
		


plot(fund.ds[,"nodf"], samp1000.ds[,"nodf"], type="n")
abline(0,1)

CI95.2D(fund.ds$nodf[fund.dd$link.rule == "trait.narrow"],
		samp1000.ds$nodf[samp1000.dd$link.rule == "trait.narrow"],
		col=hsv(0.7, alpha=0.5))
		
CI95.2D(fund.ds$nodf[fund.dd$link.rule == "trait.narrow"],
		eco.ds$nodf[eco.dd$link.rule == "trait.narrow"],
		col=hsv(0.8, alpha=0.5))
		
CI95.2D(fund.ds$nodf[fund.dd$link.rule == "trait.narrow"],
		samp100.ds$nodf[samp100.dd$link.rule == "trait.narrow"],
		col=hsv(0.9, alpha=0.5))
		








