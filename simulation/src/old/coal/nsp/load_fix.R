load.fix <- function(files, paths, nrule, nmat, ntopo){

  load.name <- function(path, file){
    load(file.path(path,file))
    return(res)
  } 
  files.load <- list.files(files)
  dat <- lapply(files.load, load.name, path= paths)
  nsim <- length(dat)*length(dat[[2]][[2]])
  parms <- lapply(dat, function(X){X[[1]]}) ## first is the params
  data <- lapply(dat, function(X){X[[2]]}) ## second element is the data
  data <- lapply(data, function(X){do.call(rbind, X)})
  data <- do.call(rbind, data)
  
  link.rule <- rownames(data)
  rownames(data) <- NULL
  
  topo <- rep(c(rep("same.same", nrule*nmat),
                rep("same.diff", nrule*nmat),
                rep("diff.diff", nrule*nmat),
                rep("diff.same", nrule*nmat)), nsim)
  
  mats <- rep(c("evo", "eco", "eco2"), nsim*nrule*ntopo)
  
  mean.abund <- sapply(parms, function (x) x$mean.abund)
  mean.abund <- rep(mean.abund, each= ntopo*nrule*nmat)
  
  sd.abund <- sapply(parms, function (x) x$sd.abund)
  sd.abund <- rep(sd.abund, each= ntopo*nrule*nmat)
  
  w.evo <- sapply(parms, function (x) x$w.evo)
  w.evo <- rep(w.evo, each= ntopo*nrule*nmat)

  nsp <- sapply(parms, function (x) x$sp)
  nsp <- rep(nsp, each= ntopo*nrule*nmat)
  
  cats <- cbind(link.rule, topo, mats)
  data <- cbind(data, mean.abund, sd.abund, w.evo, nsp)

  cats[,"link.rule"] <- as.character(cats[,"link.rule"])
  cats[,"link.rule"][grep("neutral", cats[,"link.rule"])] <- "neutral"
  cats[,"link.rule"][grep("barrior", cats[,"link.rule"])] <- "barrior"
  cats[,"link.rule"][grep("matching", cats[,"link.rule"])] <-
  "matching"
  
  return(list(data, cats))
}
