set.seed(3)
cases <- list(mu=seq(from= 0.01, to= 1.0, length.out= 10),
                     lambda=seq(from= 0.01, to= 1.0, length.out= 10))
                     
for(i in 1:length(cases[[1]])){
	print(c(cases[[1]][i],cases[[2]][i]))
	tree <- sim.bd.taxa(15, 1, lambda = cases[[1]][i], mu=cases[[1]][i], complete=FALSE)
	plot(tree[[1]][[1]]
	print(tree[[2]]))
}      