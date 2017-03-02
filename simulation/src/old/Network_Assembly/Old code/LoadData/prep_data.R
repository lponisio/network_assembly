prep.data <- function(data, data.cat, ages, nspecies){
	
	abrev.data <- data[data.cat$age == ages & sim.data.cat$npol == nspecies,]

	abrev.cat <- data.cat[data.cat$age == ages & data.cat$npol == nspecies,]

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


	return(list(fund.ss = fund.ss, fund.sd= fund.sd, fund.ds= fund.ds, fund.dd=fund.dd,
				eco.ss = eco.ss, eco.sd= eco.sd, eco.ds= eco.ds, eco.dd=eco.dd,
				samp1000.ss = samp1000.ss, samp1000.sd= samp1000.sd, samp1000.ds= samp1000.ds, samp1000.dd=samp1000.dd,
				samp100.ss = samp100.ss, samp100.sd= samp100.sd, samp100.ds= samp100.ds, samp100.dd=samp100.dd))
}