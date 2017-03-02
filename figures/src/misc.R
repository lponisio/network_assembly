pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

add.alpha <- function(col, alpha=0.2){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
        rgb(x[1], x[2], x[3],
            alpha=alpha))
}

inv.logit <- function(a){
  exp(a)/(exp(a) + 1)
}


mean.sd <- function(x){
  m <- mean(x, na.rm=TRUE)
  s <- sd(x, na.rm=TRUE)
  ci.lb <- m - 1.96*s
  ci.ub <- m + 1.96*s
  return(c(m, s, ci.lb, ci.ub))
}
