# alpine_saxifraga
Scripts used to perform analyses for: Repeated upslope biome shifting during late-Cenozoic climate cooling in a diverse alpine plant clade

Divergence time estimation: this directory provides input files for *treePL* to enable estimation of each of the different time calibrated phylogenies that were estimated in this study. 

[divergence_time_estimation](https://github.com/pebgroup/alpine_saxifraga/tree/main/divergence_time_estimation): provides the necessary inputs to estimate the time-calibrated phylogenies that were estimated in this study. In each case the *treePL* configuration file is provided (ending in .treePL) and the input tree file (ending in .tre)

[lineage_specific_diversification_rate_estimation](https://github.com/pebgroup/alpine_saxifraga/tree/main/lineage_specific_diversification_rate_estimation): provides the necessary inputs to perform lineage specific diversification rate estimation on the time-calibrated phylogeny designated as main, and no maximum. In each case the BAMM input file, and input time-calibrated phylogeny is provided. Also provided is an [R script](https://github.com/pebgroup/alpine_saxifraga/blob/main/lineage_specific_diversification_rate_estimation/set_priors.R) used to get appropriate priors in each case, and an [R script](https://github.com/pebgroup/alpine_saxifraga/blob/main/lineage_specific_diversification_rate_estimation/bamm_analysis.R) used for analysing and plotting outputs.

CLASSE: for performing classe analyses in RevBayes<br>
diversification: for estimating diversification rates in BAMM<br>
