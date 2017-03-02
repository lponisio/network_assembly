source("/Users/laurenponisio/Documents/Network_Assembly/load_fix.R")

setwd("/Users/laurenponisio/Documents/Network_Assembly")


data1 <- load.fix(files="03_12_2012",paths="03_12_2012",nrule=8,nmat=5,ntopo=4)

data2 <- load.fix(files="03_13_2012",paths="03_13_2012",nrule=8,nmat=5,ntopo=4)

#data3 <- load.fix(files="03_14_2012",paths="03_14_2012",nrule=8,nmat=5,ntopo=4)

#data3 <- data3[-grep("Error",as.character(data3[,1])),]


data4 <- load.fix(files="03_15_2012",paths="03_15_2012",nrule=8,nmat=5,ntopo=4)

#data5 <- load.fix(files="03_16_2012",paths="03_16_2012",nrule=8,nmat=5,ntopo=4)

data6 <- load.fix(files="03_17_2012",paths="03_17_2012",nrule=8,nmat=5,ntopo=4)

sim.data <- rbind(data1,data2, data4, data6)
sim.data$mean.cor <- (sim.data$cor.pol + sim.data$cor.plant)/2

sim.data.data <- sim.data[,1:9]
sim.data.cat <- sim.data[,10:17]
sim.data.data <- as.matrix(sim.data.data)



save.image("loaded_data_prelim.RData")