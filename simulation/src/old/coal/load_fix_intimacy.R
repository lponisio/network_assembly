## load data with only evo mats
load.fix <- function(files, nrule=2, nmat=1, ntopo=4){

  load.name <- function(path, file){
    load(file.path(path,file))
    return(res)
  }
  files.load <- list.files(files)
  dat <- lapply(files.load, load.name, path=files)
  nsim <- length(dat)*length(dat[[2]][[2]])
  parms <- lapply(dat, function(X){X[[1]]}) ## first is the params
  data <- lapply(dat, function(X){X[[2]]}) ## second element is the data
  data <- lapply(data, function(X){do.call(rbind, X)})
  data <- do.call(rbind, data)

  colnames(data) <- c("NODF", "H2", "modularityG", "modularityR",
                      "modularityD",
                      "zNODF", "zH2", "zmodG", "zmodR", "zmodD",
                      "pNODF", "pH2", "pmodG", "pmodR", "pmodD",
                      "corPol","corPlant", "fill",
                      "npol", "nplant", "connectance")

  topo <- rep(c(rep("same.same", nrule*nmat),
                rep("same.diff", nrule*nmat),
                rep("diff.diff", nrule*nmat),
                rep("diff.same", nrule*nmat)), nsim)
  
  range.size <- sapply(parms, function (x) x$range.size)
  range.size <- rep(range.size, each= ntopo*nrule*nmat)

  link <- do.call(rbind, strsplit(rownames(data), '\\.'))
  colnames(link) <- c("link.rule")
  mats <- rep("evo", nrow(data))
  
  rownames(data) <- NULL
  cats <- cbind(link, mats, topo)
  data <- cbind(data, range.size)
  return(list(data, cats))
}
