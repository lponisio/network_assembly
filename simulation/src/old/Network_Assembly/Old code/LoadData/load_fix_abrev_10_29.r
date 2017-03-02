
load.name <- function(path, file){
  load(file.path(path,file))
  return(res)
}

load.fix <- function(files, paths, nrule, nmat, ntopo) {

  ## **************************************************
  ## create and repair data
  files_load <- list.files(files)
  dat <- lapply(files_load, load.name, path=paths)

  data <- lapply(dat, function(X){X[[2]]}) #take the first element, the data

  f <- function(dd) {
    drop.broken <- function(x)
      x[-which(sapply(x, length)==1)]
    drop.broken(dd)
  }
  data <- lapply(data, f)
  nreps <- sapply(data, length)
  
  ## **************************************************

  ## **************************************************
  ## create parameter matrix
  parms <- t(sapply(dat, function(X){X[[1]]})) # take the second element,
                                        # the parameters 

  ## rbind all of the elements of the list 
  data <- do.call(rbind, lapply(data, do.call, what=rbind))

  data.cats <- as.data.frame(rownames(data))
  rownames(data) <- NULL
  data <- as.data.frame(data)
  data.cats <- cbind(data.cats, rep(c(rep("same.same", nrule*nmat),
                                      rep("same.diff", nrule*nmat),
                                      rep("diff.diff", nrule*nmat),
                                      rep("diff.same", nrule*nmat)),
                                    sum(nreps)))

  data.cats <- cbind(data.cats, rep(c("fund", "eco"),
                                    sum(nreps)*nrule*ntopo))

  prms.rep <- parms[rep(1:length(nreps), nreps*nrule*nmat*ntopo),
                    c("mu", "lambda", "age")]
  data.cats <- cbind(data.cats, prms.rep)
  colnames(data.cats) <- c("link.rule", "topo", "mat", "mu", "lambda",
                           "age")
  ## **************************************************

  data.cats[,"link.rule"] <- as.character(data.cats[,"link.rule"])
  data.cats[,"link.rule"][grep("neutral", data.cats[,"link.rule"])] <- "neutral"

  data.cats[,"link.rule"][grep("trait.narrow", data.cats[,"link.rule"])] <- "trait.narrow"

  data.cats[,"link.rule"][data.cats[,"link.rule"] == "trait.bar1"] <- "trait.bar"
  data.cats[,"link.rule"][data.cats[,"link.rule"] == "trait.bar2"] <- "trait.bar"

  data.cats[,"link.rule"] <- as.factor(data.cats[,"link.rule"])

  return(list(dats=data,cats=data.cats))
}
