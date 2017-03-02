setwd("/Users/laurenponisio/Documents/Network_Assembly")
source("Functions/load_fix_abrev_10_29.r")
source("Functions/prep_plot.R")

prelim.res <- load.fix(files="10_18_2012",paths="10_18_2012",nrule=3,nmat=2,ntopo=4)

prep.plot(prelim.res, "nodf")

prep.plot.by.diff(prelim.res, "nodf", legend.loc="topright", place.legend="second")

prep.plot.by.diff(prelim.res, "mod", legend.loc="topright", place.legend="first")
