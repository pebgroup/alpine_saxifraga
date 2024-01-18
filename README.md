# Repeated upslope biome shifting during late-Cenozoic climate cooling in a diverse alpine plant clade     

## [gene trees, alignments, and raw distribution data](https://zenodo.org/records/10530003)
Follow [this link](https://zenodo.org/records/10530003) to the zenodo repository.

## [molecular branch length estimation](https://github.com/pebgroup/alpine_saxifraga/tree/main/molecular_branch_length_estimation)
This contains scripts for branch-wise analysis and gene shopping. Needs a reference [species tree](https://github.com/pebgroup/alpine_saxifraga/blob/main/molecular_branch_length_estimation/species_tree_rooted.tre) and gene trees. Gene trees are available on the related Zenodo repository. The branch_wise script is an early version of the script that was used for this study. A far faster version will be made available soon. 

## [divergence time estimation](https://github.com/pebgroup/alpine_saxifraga/tree/main/divergence_time_estimation)
This provides the necessary inputs to estimate the time-calibrated phylogenies that were estimated in this study. In each case the *treePL* configuration file is provided (ending in .treePL) and the input tree file (ending in .tre)

## [lineage specific diversification rate estimation](https://github.com/pebgroup/alpine_saxifraga/tree/main/lineage_specific_diversification_rate_estimation)
This provides the necessary inputs to perform lineage specific diversification rate estimation on the time-calibrated phylogeny designated as main, and no maximum. In each case the BAMM input file, and input time-calibrated phylogeny is provided. Also provided is an [R script](https://github.com/pebgroup/alpine_saxifraga/blob/main/lineage_specific_diversification_rate_estimation/set_priors.R) used to get appropriate priors in each case, and an [R script](https://github.com/pebgroup/alpine_saxifraga/blob/main/lineage_specific_diversification_rate_estimation/bamm_analysis.R) used for analysing and plotting outputs.

## [ClaSSE](https://github.com/pebgroup/alpine_saxifraga/tree/main/ClaSSE) 
This contains all the scripts needed to perform ClaSSE analyses that were undertaken as part of this study. [main](https://github.com/pebgroup/alpine_saxifraga/tree/main/ClaSSE/main) and [no_maximum](https://github.com/pebgroup/alpine_saxifraga/tree/main/ClaSSE/no_maximum) contain scripts for analysis on the main and no maximum trees respectively. [final_tree_and_data_frame_processing.r](https://github.com/pebgroup/alpine_saxifraga/blob/main/ClaSSE/final_tree_and_data_frame_processing.r) prunes the species tree and compiles the [biome_and_region_preference](https://github.com/pebgroup/alpine_saxifraga/blob/main/ClaSSE/biome_region.csv) into a matrix for analysis in the ClaSSE model. Note that line 118 in [final_tree_and_data_frame_processing.r](https://github.com/pebgroup/alpine_saxifraga/blob/main/ClaSSE/final_tree_and_data_frame_processing.r) needs adjusting depending on which output tree from treePL you are working on. [biome_and_region_preference](https://github.com/pebgroup/alpine_saxifraga/blob/main/ClaSSE/biome_region.csv) is equivalent to Table S5 in the manuscript.  

## [Generate Saxifraga distribution map in Figure 1](https://github.com/pebgroup/alpine_saxifraga/tree/main/Generate_Figure_1)
Scripts to generate the distribution map in figure 1. Source data for the distribution points is [here](https://zenodo.org/records/8408326)

## [Generate Plots in Figure 2](https://github.com/pebgroup/alpine_saxifraga/tree/main/Generate_Figure_2)
Scripts for plotting output from BAMM analysis (Fig. 2a and b) and biome specific speciation rates from CLASSE (Fig. 2c). Log files from BAMM and Classe [here](https://zenodo.org/records/8408326). 

## [Generate Plots in Figure 3](https://github.com/pebgroup/alpine_saxifraga/tree/main/Generate_Figure_3)
Scripts for plotting rates of upslode biome shifts and inter regional mountain hopping as displayed in Figure 3. Simmaps are [here](https://zenodo.org/records/8408326)

