## load data with only evo mats
load.fix <- function(files, nrule=prms$nrule, nmat=prms$nmat,
                     ntopo=4, bd){

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
  colnames(data) <- c('NODF',
                      'modularityG', 'modularityR','modularityD',
                      'zNODF', 'zmodG', 'zmodR', 'zmodD',
                      'pNODF', 'pmodG', 'pmodR', 'pmodD',
                      'corPol','corPlant',
                      'meanDisPol', 'meanDisPlant',
                      'fill',
                      'npol', 'nplant',
                      'meanIntPol', 'meanIntPlant',
                      'connectance')

  rconnectance <- round(data[,'connectance'], 1)
  data <- cbind(data, rconnectance)

  ## pp ratio by metrics
  ratio <- round(log(data[,'npol']) - log(data[,'nplant']), 1)
  data <- cbind(data, ratio)

  ## mean proportion
  rdegree <- round(data[,'meanIntPol'], 1)
  data <- cbind(data, rdegree)

  ## mean similarity
  rdis <- round(data[,'meanDisPol'], 1)
  data <- cbind(data, rdis)

  ## log fill
  logfill <- log(data[,'fill'] + 10^-10)
  data <- cbind(data, logfill)

  ## nsp by metrics
  nsp <- round(((data[,'npol'] + data[,'nplant'])/5))*5
  data <- cbind(data, nsp)

  topo <- rep(c(rep('same.same', nrule*nmat),
                rep('same.diff', nrule*nmat),
                rep('diff.diff', nrule*nmat),
                rep('diff.same', nrule*nmat)), nsim)
  
  range.size <- sapply(parms, function (x) x$range.size)
  range.size <- rep(range.size, each= ntopo*nrule*nmat)

  link <- do.call(rbind, strsplit(rownames(data), '\\.'))
  colnames(link) <- c('link.rule')
  mats <- rep('evo', nrow(data))
  
  rownames(data) <- NULL
  cats <- cbind(link, mats, topo)
  if(bd == 'bd'){
    mu <- sapply(parms, function (x) x$mu)
    mu <- rep(mu, each= ntopo*nrule*nmat)
    
    lambda <- sapply(parms, function (x) x$lambda)
    lambda <- rep(lambda, each= ntopo*nrule*nmat)

    age <- sapply(parms, function (x) x$age)
    age <- rep(age, each= ntopo*nrule*nmat)

    data <- cbind(data,
                  range.size,
                  mu,
                  lambda,
                  age)

  } else{
    data <- cbind(data,
                  range.size)
  }
  return(list(data, cats))
}
