
library(BAMMtools)
library(coda)

###GET_EVENT_DATA_AND_BAMM_DATA_OBJECT

tree <- read.tree(...)
event_data <- getEventData(tree, eventdata=..., burnin=0.1, nsamples=2000)

###ASSESS_CONVERGENCE

mcmcout <- read.csv(..., header=T) ###USES_MCMC_OUTPUT_FILE
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout)) ###burnin_for_mcmcout
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

ess_n_shifts <- effectiveSize(postburn$N_shifts)
ess_likelihood <- effectiveSize(postburn$logLik)

###GET_POST_PROBS_OF_DIFFERENT_NUMBERS_OF_SHIFTS

post_probs <- table(postburn$N_shifts) / nrow(postburn)

###SAME

shift_probs <- summary(event_data)

###CALCULATE_BAYES_FACTORS_FOR_DIFFERENT_MODELS

bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1)

###GET_MEAN_PHYLORATE_PLOT

plot.bammdata(event_data, lwd=2, legend=T)

###GET_CREDIBLE_SHIFT_SET

css <- credibleShiftSet(event_data, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95, legend=T)

###GET_BEST_TREE

best <- getBestShiftConfiguration(event_data, expectedNumberOfShifts=1)

###FINAL_PLOTTING

plot(best, labels = TRUE, cex= 0.3, lwd=2, legend=T)
plot.credibleshiftset(css)
