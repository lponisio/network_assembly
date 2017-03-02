rm(list=ls())
setwd("/Users/laurenponisio/Documents/Network_Assembly")
files <- "10_18_2012"
paths <- "10_18_2012"
nrule <- 3
nmat <- 2
ntopo <- 4
nrep <- 100

load.name <- function(path, file){
	load(file.path(path,file))
	return(res)
}




load.fix <- function(files, paths, nrule, nmat, ntopo, nrep){

	files_load <- list.files(files)
	dat <- lapply(files_load, load.name, path=paths)
	
	data <- lapply(dat, function(X){X[[2]]}) #take the first element, the data
	
	f <- function(dd) {
		drop.broken <- function(x)
			x[-which(sapply(x, length)==1)]
		drop.broken(dd)
	}
	data.fixed <- lapply(data, f)
	nreps <- sapply(data.fixed, length)

	ntree <- length(dat)*length(dat[[2]][[2]]) ## number of total rows

	parms <- lapply(dat, function(X){X[[1]]}) # take the second element, the parameters
	
	data <- lapply(data, function(X){do.call(rbind, X)}) # rbind all of the elements of the list

	data <- do.call(rbind, data) # rbind again
	
	data.cats <- as.data.frame(rownames(data))
	rownames(data) <- NULL
	data <- as.data.frame(data)
	data$cats <- rep(1:100, nreps)

2*3*4*sum(nreps)


	data.cats <-cbind(data.cats, rep(c(rep("same.same", nrule*nmat), rep("same.diff", nrule*nmat), rep("diff.diff", nrule*nmat), rep("diff.same", nrule*nmat)), ntree)[1:nrow(data)])

	data.cats<-cbind(data.cats, rep(c("fund", "eco"), ntree*nrule*ntopo)[1:nrow(data.cats)])

	mus <- 0
	lambdas <- 0
	ages <- 0
	
	for(i in 1:length(parms)){
		 mus <- c(mus, rep(parms[[i]]$mu, nrep*nrule*nmat))
		lambdas <- c(lambdas, rep(parms[[i]]$lambda, nrep*nrule*nmat))
		ages <- c(ages, rep(parms[[i]]$age, nrep*nrule*nmat))
	}
	mus <- mus[-1]
	lambdas <- lambdas[-1]
	ages <- ages[-1]
	
	data.cats <- cbind(data.cats, mus, lambdas, ages)
	colnames(data.cats) <- c("link.rule", "topo", "mat", "mu", "lambda", "age")

	data.cats[,"link.rule"] <- as.character(data.cats[,"link.rule"])
	data.cats[,"link.rule"][grep("neutral", data.cats[,"link.rule"])] <- "neutral"
	
	data.cats[,"link.rule"][grep("trait.narrow", data.cats[,"link.rule"])] <- "trait.narrow"

	data.cats[,"link.rule"][data.cats[,"link.rule"] == "trait.bar1"] <- "trait.bar"
	data.cats[,"link.rule"][data.cats[,"link.rule"] == "trait.bar2"] <- "trait.bar"

	data.cats[,"link.rule"] <- as.factor(data.cats[,"link.rule"])
	
	return(list(dats=data,cats=data.cats))
}