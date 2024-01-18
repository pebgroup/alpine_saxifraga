library(ggplot2)
library(BAMMtools)
library(coda)
library(phytools)
library(devtools)
library(RRphylo)
library(phangorn)

###GET_EVENT_DATA_AND_BAMM_DATA_OBJECT

tree <- read.tree("analysis_tree.tre")
event_data <- getEventData(tree, eventdata="sax_div_data.txt", burnin=0.1, nsamples=2000)

###GET_CREDIBLE_SHIFT_SET

css <- credibleShiftSet(event_data, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95, legend=T)
plot.credibleshiftset(css, spex="netdiv", legend=T)

###GET_BEST_TREE

best <- getBestShiftConfiguration(event_data, expectedNumberOfShifts=1)
ggsave("bamm_plot.pdf",plot(best, labels = TRUE, cex=0.125, lwd=0.7, legend=T) + addBAMMshifts(best, cex=1.5) + add.scale.bar(),limitsize=FALSE)

###PLOT RATE THROUGH TIME

plotRateThroughTime(event_data, intervals = seq(from = 0.05, to = 0.95, by = 0.01), useMedian=FALSE, ratetype='auto')