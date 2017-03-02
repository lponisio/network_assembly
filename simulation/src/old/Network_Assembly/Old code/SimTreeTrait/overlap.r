rm(list=ls())
setwd("~/Desktop")
data <- read.csv("watertemp.csv", as.is=TRUE)

species <- data$Species
temps <- data[c("Mintemp", "Maxtemp")]
temps[temps==100] <- NA

f <- function(range1, range2) {
  if(any(is.na(c(range1, range2)))) return(NA)
  min1 <- min(range1)
  max1 <- max(range1)
  min2 <- min(range2)
  max2 <- max(range2)
  if(max1<=min2 | max2<=min1)
    return(0)
  if((max1>=max2 & min1<=min2) | (max2>=max1 & min2<=min1))
    return(1)
  if(max1>=max2 & min1>=min2)
    return((max2-min1)/min(max1-min1, max2-min2))
  if(max2>=max1 & min2>=min1)
    return((max1-min2)/min(max1-min1, max2-min2))
}

temp.ind <- expand.grid(1:nrow(data), 1:nrow(data))

m <- mapply(function(a, b) f(temps[a,], temps[b,]),
            a=temp.ind[,1], b=temp.ind[,2])
m <- matrix(m, ncol=nrow(data), nrow=nrow(data))
rownames(m) <- species
colnames(m) <- species
