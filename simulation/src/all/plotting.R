library(RColorBrewer)
options(warn=-1)
if(tree == "bd") {
  range.intint.quan.spec <- c(0.1, 0.9)
  range.intint.qual.spec <- c(0.1, 0.9)

  range.intint.qual.ol <- c(0.5, 0.9)
  range.intint.quan.ol <- c(0.5, 0.9)

} else {
  range.intint.quan.spec <- c(0, 0.9)
  range.intint.qual.spec <- c(0.1, 0.9)

  range.intint.qual.ol <- c(0.5, 0.9)
  range.intint.quan.ol <- c(0.5, 0.9)

}

path.fig1 <- file.path(fig.dir, sprintf('intimacy/%s', type))
path.fig2 <- file.path(fig.dir, sprintf('NODFbyMod/%s', type))
path.fig3 <- file.path(fig.dir, sprintf('SingleMetrics/%s', type))

## interaction intimacy
## specialization diff
intInt.diff.plot3(simres=dats,
                  metric1='zmodD',
                  metric2='zNODF',
                  metric3='corPol',
                  met1lab='Relative \n Modularity',
                  met2lab='Relative \n Nestedness',
                  met3lab='Phylogenetic \n Interaction Signal',
                  xmetric='rdegree',
                  xlabel='Interaction Niche Breath',
                  path.fig=path.fig1,
                  range.intint.quan=range.intint.quan.spec,
                  range.intint.qual=range.intint.qual.spec)

## intInt.diff.image(simres=dats,
##                   metric1='zNODF',
##                   xmetric1='rdegree',
##                   xmetric2='rdis',
##                   path.fig=path.fig1)

## intInt.diff.image(simres=dats,
##                   metric1='zmodD',
##                   xmetric1='rdegree',
##                   xmetric2='rdis',
##                   path.fig=path.fig1)

## intInt.diff.image(simres=dats,
##                   metric1='corPol',
##                   xmetric1='rdegree',
##                   xmetric2='rdis',
##                   path.fig=path.fig1)



## overlap diff
intInt.diff.plot3(simres=dats,
                  metric1='zmodD',
                  metric2='zNODF',
                  metric3='corPol',
                  met1lab='Relative \n Modularity',
                  met2lab='Relative \n Nestedness',
                  met3lab='Phylogenetic \n Interaction Signal',
                  xmetric='rdis',
                  xlabel='Partner overlap',
                  path.fig=path.fig1,
                  range.intint.qual=range.intint.qual.ol,
                  range.intint.quan=range.intint.quan.ol,
                  right.lab="High",
                  left.lab="Low")

## specialization
intInt.plot3(simres=dats,
             metric1='zmodD',
             metric2='zNODF',
             metric3='corPol',
             met1lab='Relative \n Modularity',
             met2lab='Relative \n Nestedness',
             met3lab='Phylogenetic \n Interaction Signal',
             xmetric='rdegree',
             xlabel='Interaction Niche Breath',
             path.dir=path.fig1,
             rev.x=FALSE,
             leg.panel='left',
             range.intint=c(0, 0.9))

## overlap
intInt.plot3(simres=dats,
             metric1='zmodD',
             metric2='zNODF',
             metric3='corPol',
             met1lab='Relative \n Modularity',
             met2lab='Relative \n Nestedness',
             met3lab='Phylogenetic \n Interaction Signal',
             xmetric='rdis',
             xlabel='Partner overlap',
             leg.loc='bottomleft',
             leg.panel='left',
             path.dir=path.fig1,
             range.intint=c(0.5, 0.9),
             right.lab="High",
             left.lab="Low")

## overlap and specialization by intint
intInt.plot2(simres=dats,
             metric1='rdegree',
             metric2='rdis',
             met1lab='Interaction \n Niche Breath',
             met2lab='Partner overlap',
             xmetric='range.size',
             xlabel='Trait range',
             leg.panel="right",
             leg.loc='bottomleft',
             path.dir=path.fig1,
             left.lab="Narrow",
             right.lab="Wide")


## nestedness by modularity bar plot

mets3.bar(simres = dats,
          path.fig=path.fig2,
          metric1='zmodD',
          metric2='zmodR',
          metric3='zmodG',
          column='range.size',
          subset=FALSE,
          met1lab='Edge Betweeness',
          met2lab='Random Walk',
          met3lab='Greedy',
          adj.lab=1)

mets3.bar(simres = dats,
          path.fig=path.fig2,
          metric1= 'zmodD',
          metric2='zNODF',
          metric3= 'corPol',
          column='range.size',
          subset=FALSE,
          met3lab='Phylogenetic \n Interaction Signal',
          adj.lab=0.1)


diff.bar3(simres=dats,
          metric1='zmodD',
          metric2='zNODF',
          metric3='corPol',
          met1lab='Relative \n Modularity',
          met2lab='Relative \n Nestedness',
          met3lab='Phylogenetic \n Interaction Signal',
          path.fig=path.fig2,
          adj.lab=0.8)



## connectance by metrics
## intInt.plot2(simres=dats,
##              metric1='zmodD',
##              metric2='zNODF',
##              met1lab='Modularity',
##              met2lab=' Nestedness',
##              xmetric='connectance',
##              xlabel='Connectance',
##              leg.loc='bottomleft',
##              path.dir=path.fig3)

## intInt.plot2(simres=dats,
##              metric1='zmodD',
##              metric2='zNODF',
##              met1lab='Modularity',
##              met2lab=' Nestedness',
##              xmetric='nsp',
##              xlabel='Number of Species',
##              leg.loc='bottomright',
##              path.dir=path.fig3)



## intInt.plot2(simres=dats,
##              metric1='zmodD',
##              metric2='zNODF',
##              met1lab='Modularity',
##              met2lab=' Nestedness',
##              xmetric='ratio',
##              xlabel='Species ratio (log)',
##              leg.loc='bottomright',
##              path.dir=path.fig3,
##              drop.same.same=TRUE)


