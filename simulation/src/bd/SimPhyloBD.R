sim.phylo <- function(mu,
                      lambda,
                      age,
                      nspecies) {
    ## takes parameters for a birth death tree and generates two trees

  trees <- sim.bd.taxa.age(n=nspecies,
                            numbsim=2,
                            lambda=lambda,
                            mu= mu,
                            frac= 1,
                            age= age,
                            mrca= FALSE)

  return(list(tree.1=trees[[1]], tree.2=trees[[2]]))
}

