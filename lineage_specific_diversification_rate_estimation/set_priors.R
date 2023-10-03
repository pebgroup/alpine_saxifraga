library(phytools)
library(BAMMtools)

tree <- read.tree("...")
setBAMMpriors(tree, total.taxa = 557, traits = NULL, outfile = "myPriors.txt", Nmax = 1000, suppressWarning = FALSE)