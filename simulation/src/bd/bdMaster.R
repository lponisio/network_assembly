run.sim <- function(prms, i, save.path) {
    ## function to generate network communities using bd trees
    prms <- prms.f(prms)
    prms$mu <- cases[i,'mu']
    prms$lambda <- cases[i,'lambda']
    prms$range.size <- cases[i, 'range.size']

    all.trees <- sim.phylo(mu=prms$mu,
                           lambda=prms$lambda,
                           nspecies=prms$sp,
                           age=prms$age)

    f <- function(x) {
        master.fun(prms,
                   tree1 = all.trees[[1]],
                   tree2 = all.trees[[2]],
                   nnul = prms$nulls)
    }
    print(i)
    sim.links <- lapply(1:prms$ntree, f)
    res <- list(prms= prms, sim.links= sim.links)
    save(res, file= file.path(save.path, sprintf('%d.RData', i)))
}

