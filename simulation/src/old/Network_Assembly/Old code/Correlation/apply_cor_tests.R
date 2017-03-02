apply.cor.test <- function(data, metrics){
	cor.tests(
	xdatas=list(ss=data$fund.ss, ds=data$fund.ds, sd=data$fund.sd,dd=data$fund.dd), 
	ydata1s=list(data$eco.ss,data$eco.ds,data$eco.sd,data$eco.dd), 
	ydata2s=list(data$samp1000.ss,data$samp1000.ds,data$samp1000.sd,data$samp1000.dd), 
	ydata3s=list(data$samp100.ss,data$samp100.ds,data$samp100.sd,data$samp100.dd), trait=c("trait.narrow", "trait.bar", "neutral"), metric=metrics)

}