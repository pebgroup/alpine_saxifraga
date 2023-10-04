library(ape)
library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)

overall_species_tree <- read.tree("species_tree_rooted.tre")
error_bootstrap_threshold <- 75
error_bootstrap_proportion_threshold <- 0.1
min_coverage <- 2

gene_tree_name <- vector("list", 0)
gene_tree_index <- vector(mode="numeric", length=0)

wrong_topology_support_proportion <- vector(mode="numeric", length=0)
number_of_taxa_for_loci <- vector(mode="numeric", length=0)
proportion_for_loci <- vector(mode="numeric", length=0)

######################
###DEFINE_OUTGROUPS###
######################

outgroup_taxa <- vector("list", 12)
outgroup_taxa[[1]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Medusandra_richardsiana_P05_WA02", "Soyauxia_talbotii_P01_WC04"), "node"))$tip.label
outgroup_taxa[[2]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Paeonia_morisii_P03_WH11", "Distyliopsis_laurifolia_P01_WA04"), "node"))$tip.label
outgroup_taxa[[3]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Penthorum_chinense_P01_WC06", "Sedum_filipes_P01_WA08")))$tip.label
outgroup_taxa[[4]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Pterostemon_mexicanus_P01_WC01", "Itea_chinensis_S371")))$tip.label
outgroup_taxa[[5]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Ribes_griffithii_P06_WB04", "Ribes_acuminatum_P06_WE06")))$tip.label
outgroup_taxa[[6]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Leptarrhena_pyrolifolia_P01_WC10", "Micranthes_punctata_S251")))$tip.label
outgroup_taxa[[7]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Cascadia_nuttallii_P06_WB12", "Saxifragodes_albowiana_S375")))$tip.label
outgroup_taxa[[8]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Saxifraga_mertensiana_P07_WF04", "Saxifraga_imparilis_P03_WC02")))$tip.label
outgroup_taxa[[9]] <- "Saxifragella_bicuspidata_P01_WD08"
outgroup_taxa[[10]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Saxifraga_lactea_S260", "Saxifraga_taylorii_P05_WE05")))$tip.label
outgroup_taxa[[11]] <- "Saxifraga_josephii_S345"
outgroup_taxa[[12]] <- extract.clade(overall_species_tree, findMRCA(overall_species_tree, c("Saxifraga_eschscholtzii_S277", "Saxifraga_consanguinea_P06_WB01")))$tip.label

unrootable <- vector(mode="numeric", length=0)

############
###REROOT###
############

`%notin%` <- Negate(`%in%`)
n_gene_trees <- 329
gene_trees <- vector("list", n_gene_trees)
species_trees <- vector("list", n_gene_trees)
setwd("RAXML_OUTPUT/")
for (a in 1:329){
print(paste("rooting locus", a, sep=""))	 
gene_trees[[a]] <- read.tree(list.files(pattern = "bipartitions.Locus")[[a]])
setwd("..")
species_trees[[a]] <- read.tree("species_tree.tre")
setwd("RAXML_OUTPUT/")
species_trees[[a]] <- drop.tip(species_trees[[a]], which(species_trees[[a]]$tip.label %notin% gene_trees[[a]]$tip.label))  
keep_searching <- "YES"
for (i in 1:length(outgroup_taxa)){
if (keep_searching =="YES"){
if (length(which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)) > 1){
gene_trees[[a]] <- reroot(gene_trees[[a]], node=findMRCA(gene_trees[[a]], outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)], "node"), 0.01)
species_trees[[a]] <- reroot(species_trees[[a]], node=findMRCA(species_trees[[a]], outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% species_trees[[a]]$tip.label)], "node"), 0.01)
keep_searching <- "no"
}
else if (length(which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)) == 1){
gene_trees[[a]] <- reroot(gene_trees[[a]], which(gene_trees[[a]]$tip.label == outgroup_taxa[[i]][[which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)]]), 0.01) 
species_trees[[a]] <- reroot(species_trees[[a]], which(species_trees[[a]]$tip.label == outgroup_taxa[[i]][[which(outgroup_taxa[[i]] %in% species_trees[[a]]$tip.label)]]), 0.01)
keep_searching <- "no"
}
}
}
if (keep_searching == "YES"){
unrootable <- append(unrootable, a) 
keep_searching <- "no"
}
}

##########################
###PERFORM_KEY_ANALYSES###
##########################

for (a in 1:329){
print(paste("performing gene tree analysis", a, sep="_"))
wrong_topology_support_counter <- vector(mode="numeric", length=0)
species_tree_clades <- vector("list", (length(species_trees[[a]]$tip.label) - 3))
gene_tree_clades <- vector("list", (length(gene_trees[[a]]$tip.label) - 3))
for (j in 1:length(gene_tree_clades)){
species_tree_clades[[j]] <- extract.clade(species_trees[[a]], length(species_trees[[a]]$tip.label) + 1 + j)
gene_tree_clades[[j]] <- extract.clade(gene_trees[[a]], length(species_trees[[a]]$tip.label) + 1 + j)
}

#############
#############

removal <- vector(mode="numeric", length=0)
for (j in 1:length(gene_tree_clades)){
if (as.numeric(gene_tree_clades[[j]]$node.label[[1]]) %in% seq(1, 100, 1)){
if (as.numeric(gene_tree_clades[[j]]$node.label[[1]]) < 0.75){
removal <- append(removal, j)
}
}
}

if (length(removal > 0)){
gene_tree_clades <- gene_tree_clades[-removal]
}

############
############

gene_tree_clades_in <- vector("list", length(gene_tree_clades))

for (j in 1:length(gene_tree_clades)){
for (i in 1:length(species_tree_clades)){
if (setequal(gene_tree_clades[[j]]$tip.label, species_tree_clades[[i]]$tip.label) == TRUE){
gene_tree_clades_in[[j]] <- 1
}
}
}
for (j in 1:length(gene_tree_clades_in)){
if (length(gene_tree_clades_in[[j]]) == 0){
if (as.numeric(gene_tree_clades[[j]]$node.label[[1]]) %in% seq(1, 100, 1)){
if (as.numeric(gene_tree_clades[[j]]$node.label[[1]]) > error_bootstrap_threshold){
wrong_topology_support_counter <- append(wrong_topology_support_counter, 1)
}
}
}
}
wrong_topology_support_proportion <- append(wrong_topology_support_proportion, length(wrong_topology_support_counter)/(length(gene_tree_clades)-1))
gene_tree_clades_in <- unlist(gene_tree_clades_in)
proportion_for_loci <- append(proportion_for_loci, (length(gene_tree_clades_in)-1)/(length(gene_tree_clades)-1))
number_of_taxa_for_loci <- append(number_of_taxa_for_loci, length(gene_trees[[a]]$tip.label))
gene_tree_name[[length(gene_tree_name)+1]] <- list.files(pattern = "bipartitions.Locus")[[a]]
gene_tree_index <- append(gene_tree_index, a)
}

##############################
###POST_ANALYSIS_PROCESSING###
##############################

search_data_frame <- data.frame(proportion_for_loci, unlist(gene_tree_name), gene_tree_index, wrong_topology_support_proportion, number_of_taxa_for_loci)
search_data_frame <- search_data_frame[order(-search_data_frame[,1]),]

removal <- vector(mode="numeric", length=0)
for (i in 1:nrow(search_data_frame)){
if (search_data_frame[,4][[i]] >  error_bootstrap_proportion_threshold){
removal <- append(removal, i)
}
}

if (length(removal) > 0){
search_data_frame <- search_data_frame[-c(removal),]
}

removal <- vector(mode="numeric", length=0)
for (i in 1:nrow(search_data_frame)){
if (search_data_frame[,5][[i]] <  0.75*length(overall_species_tree$tip.label)){
removal <- append(removal, i)
}
}

if (length(removal) > 0){
search_data_frame <- search_data_frame[-c(removal),]
}

subset_check <- vector("list", length(overall_species_tree$tip.label))
continue <- "continue"
for (j in 1:nrow(search_data_frame)){
if (continue == "continue"){
for (i in 1:length(gene_trees[[search_data_frame[,3][[j]]]]$tip.label)){
subset_check[[which(overall_species_tree$tip.label == gene_trees[[search_data_frame[,3][[j]]]]$tip.label[[i]])]] <- append(subset_check[[which(overall_species_tree$tip.label == gene_trees[[search_data_frame[,3][[j]]]]$tip.label[[i]])]], 1)
}
subset_check_lengths <- vector("list", length(subset_check))
for (i in 1:length(subset_check_lengths)){
subset_check_lengths[[i]] <- length(subset_check[[i]])
}
subset_check_lengths <- unlist(subset_check_lengths)
if (min(unlist(subset_check_lengths)) >= min_coverage){
continue <- "stop"
print(paste("use loci up to row", j, sep=""))
}
if (j == nrow(search_data_frame)){
continue <- "stop"
print("scanned all loci and not full coverage")
}
}
}

write.table(search_data_frame, "rerooted_search_data_frame_final.tsv")

for (i in 1:length(gene_trees)){
write.tree(gene_trees[[i]], paste("re_rooted_", list.files(pattern = "bipartitions.Locus")[[i]], ".tre", sep=""))
}



