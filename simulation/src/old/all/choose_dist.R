ch.dist <- function(scores, cats){
  aics <- scores[apply(scores, 1, FUN = function(x)
                       all(is.finite(x))),]
  aics.cats <- cats[apply(scores, 1, FUN = function(x)
                          all(is.finite(x))),]
  
  mins <- apply(aics, 1, min)

  for(i in 1:nrow(aics)){
    for(j in 1:ncol(aics)){
      aics[i,j] <- ifelse(aics[i,j] == mins[i], 1, 0)
    }
  }

  wins <- aggregate(aics, list(topo=aics.cats[,"topo"],
                               mats=aics.cats[,"mats"],
                               link.rule=aics.cats[,"link.rule"]),
                    function(x) sum(x)/length(x))
}
