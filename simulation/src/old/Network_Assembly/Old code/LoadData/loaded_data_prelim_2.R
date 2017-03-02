
setwd("/Users/laurenponisio/Documents/Network_Assembly")

source("Functions/load_fix.R")
source("Functions/prep_data.R")
source("Functions/cor_tests.R")
source("Functions/cor_plots_v4.R")
source("Functions/CI95_2D.R")
source("Functions/apply_cor_tests.R")
source("Functions/nodf_v_mod_plot_bymatrix.R")
source("Functions/nodf_v_mod_plot.R")
source("Functions/plot_rho.R")

data1 <- load.fix(files="04_17_2012",paths="04_17_2012",nrule=8,nmat=5,ntopo=4)

data2 <- load.fix(files="04_22_2012",paths="04_22_2012",nrule=8,nmat=5,ntopo=4)

data3 <- load.fix(files="04_24_2012",paths="04_24_2012",nrule=8,nmat=5,ntopo=4)

data4 <- load.fix(files="04_26_2012",paths="04_26_2012",nrule=8,nmat=5,ntopo=4)

data5 <- load.fix(files="05_03_2012",paths="05_03_2012",nrule=8,nmat=5,ntopo=4)

sim.data <- rbind(data1, data2, data3, data4, data5)

sim.data <- data1
sim.data$mean.cor <- (sim.data$cor.pol + sim.data$cor.plant)/2

sim.data.data <- sim.data[,1:9]
sim.data.cat <- sim.data[,10:17]
sim.data.data <- as.matrix(sim.data.data)

save.image("loaded_data_prelim_3.RData")

load("loaded_data_prelim_3.RData")

dat.15.100 <- prep.data(sim.data.data, sim.data.cat, ages=100, nspecies=15)

dat.30.100 <- prep.data(sim.data.data, sim.data.cat, ages=100, nspecies=30)

dat.60.100 <- prep.data(sim.data.data, sim.data.cat, ages=100, nspecies=60)

dat.30.50 <- prep.data(sim.data.data, sim.data.cat, ages=50, nspecies=30)

dat.30.200 <- prep.data(sim.data.data, sim.data.cat, ages=200, nspecies=30)


##calculate correlation coefficients

cors.nodf.15.100 <-apply.cor.test(dat.15.100, "nodf")

cors.nodf.30.100 <-apply.cor.test(dat.30.100, "nodf")

cors.nodf.60.100 <- apply.cor.test(dat.60.100, "nodf")

cors.nodf.30.50 <- apply.cor.test(dat.30.50, "nodf")

cors.nodf.30.200 <- apply.cor.test(dat.30.200, "nodf")



cors.mod.30.100 <- apply.cor.test(dat.30.100, "mod")


cors.mod.30.50 <- cor.tests(dat.30.50, "mod")

cors.mod.60.100 <- apply.cor.test(dat.60.100, "mod")


##plot 	
cor.plots(
	ydata1s=list(ss=dat.30.50$fund.ss, ds=dat.30.50$fund.ds, sd=dat.30.50$fund.sd,dd=dat.30.50$fund.dd), 
	xdatas=list(dat.30.50$eco.ss,dat.30.50$eco.ds,dat.30.50$eco.sd,dat.30.50$eco.dd), 
	ydata2s=list(dat.30.50$samp1000.ss,dat.30.50$samp1000.ds,dat.30.50$samp1000.sd,dat.30.50$samp1000.dd), 
	ydata3s=list(dat.30.50$samp100.ss,dat.30.50$samp100.ds,dat.30.50$samp100.sd,dat.30.50$samp100.dd), trait=c("trait.narrow", "trait.bar", "neutral"))
	
quartz()
cor.plots(
	ydata1s=list(ss=dat.30.100$fund.ss, ds=dat.30.100$fund.ds, sd=dat.30.100$fund.sd,dd=dat.30.100$fund.dd), 
	xdatas=list(dat.30.100$eco.ss,dat.30.100$eco.ds,dat.30.100$eco.sd,dat.30.100$eco.dd), 
	ydata2s=list(dat.30.100$samp1000.ss,dat.30.100$samp1000.ds,dat.30.100$samp1000.sd,dat.30.100$samp1000.dd), 
	ydata3s=list(dat.30.100$samp100.ss,dat.30.100$samp100.ds,dat.30.100$samp100.sd,dat.30.100$samp100.dd), trait=c("trait.narrow", "trait.bar", "neutral"))

modVnodf.forNSF(
	funds=list(ss=dat.30.100$fund.ss, ds=dat.30.100$fund.ds, sd=dat.30.100$fund.sd,dd=dat.30.100$fund.dd), 
	ecos=list(dat.30.100$eco.ss,dat.30.100$eco.ds,dat.30.100$eco.sd,dat.30.100$eco.dd), 
	samp1000s=list(dat.30.100$samp1000.ss,dat.30.100$samp1000.ds,dat.30.100$samp1000.sd,dat.30.100$samp1000.dd), 
	samp100s=list(dat.30.100$samp100.ss,dat.30.100$samp100.ds,dat.30.100$samp100.sd,dat.30.100$samp100.dd), trait=c("trait.narrow", "trait.bar", "neutral"))
	
modVnodf.plots.bymat(
	dat.ss= list(dat.30.100$fund.ss, dat.30.100$eco.ss, dat.30.100$samp1000.ss, dat.30.100$samp100.ss),
	dat.ds= list(dat.30.100$fund.ds, dat.30.100$eco.ds, dat.30.100$samp1000.ds, dat.30.100$samp100.ds),
	dat.sd= list(dat.30.100$fund.sd, dat.30.100$eco.sd, dat.30.100$samp1000.sd, dat.30.100$samp100.sd),
	dat.dd= list(dat.30.100$fund.dd, dat.30.100$eco.dd, dat.30.100$samp1000.dd, dat.30.100$samp100.dd), trait=c("trait.narrow", "trait.bar", "neutral"))



##plot correlations
plot.rho(cors.nodf.30.100,cors.mod.30.100)
quartz()
plot.rho(cors.nodf.30.50,cors.mod.30.50)

plot.rho(cors.nodf.60.100,cors.mod.60.100)
	