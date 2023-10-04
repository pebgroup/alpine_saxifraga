library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(coda)
library(ggplot2)
library(Rmisc)

#######################
###OVERALL_OUTGROUPS###
#######################

overall_outgroups <- c("Medusandra_richardsiana_P05_WA02", "Soyauxia_talbotii_P01_WC04", "Peridiscus_lucidus_P01_WC03")

overall_species_tree <- read.tree("species_tree_rooted.tre")
n_branches <- ((2*length(overall_species_tree$tip.label))-2)
overall_species_tree_branch_estimates <- vector("list", n_branches)

overall_species_tree_clades <- vector("list", length(overall_species_tree$tip.label) - 2)
for (i in 1:length(overall_species_tree_clades)){
overall_species_tree_clades[[i]] <- extract.clade(overall_species_tree, length(overall_species_tree$tip.label)+1+i)
}

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
outgroup_names <- vector("list", n_gene_trees)
setwd("RAXML_OUTPUT/")
for (a in 1:329){
print(paste("rooting locus", a, sep=""))	 
gene_trees[[a]] <- read.tree(list.files(pattern = "bipartitions.Locus")[[a]])
setwd("..")
setwd("RAXML_OUTPUT/")
keep_searching <- "YES"
for (i in 1:length(outgroup_taxa)){
if (keep_searching =="YES"){
if (length(which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)) > 1){
gene_trees[[a]] <- root(gene_trees[[a]], node=findMRCA(gene_trees[[a]], outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)], "node"), edgelabel=TRUE)
outgroup_names[[a]] <- outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)]
keep_searching <- "no"
}
else if (length(which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)) == 1){
gene_trees[[a]] <- root(gene_trees[[a]], outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)], , edgelabel=TRUE)
outgroup_names[[a]] <- outgroup_taxa[[i]][which(outgroup_taxa[[i]] %in% gene_trees[[a]]$tip.label)]
keep_searching <- "no"
}
}
}
if (keep_searching == "YES"){
unrootable <- append(unrootable, a) 
print("unrootable")
keep_searching <- "no"
}
}


######################
###DO_CORE_ANALYSIS###
######################

species_tree_branch_estimate <- vector("list", nrow(overall_species_tree[[1]]))

for (a in 1:length(gene_trees)){
print(a)

if (overall_outgroups[[1]] %in% gene_trees[[a]]$tip.label || overall_outgroups[[2]] %in% gene_trees[[a]]$tip.label || overall_outgroups[[3]] %in% gene_trees[[a]]$tip.label){

for (b in 1:nrow(gene_trees[[a]][[1]])){

if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
gene_trees_descendant <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,2][[b]])
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
} else {
gene_trees_descendant <- gene_trees[[a]]$tip.label[[gene_trees[[a]][[1]][,2][[b]]]]
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
}

if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
species_tree_descendant <- findMRCA(overall_species_tree, gene_trees_descendant$tip.label, "node")
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
} else {
species_tree_descendant <- which(overall_species_tree$tip.label == gene_trees_descendant)
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
}

for (c in 1:nrow(overall_species_tree[[1]])){
if (overall_species_tree[[1]][,2][[c]] == species_tree_descendant & overall_species_tree[[1]][,1][[c]] == species_tree_ancestral){
species_tree_branch_estimate[[c]] <- append(species_tree_branch_estimate[[c]], gene_trees[[a]]$edge.length[[b]])
}
}

}
} else {

ignore <- which(gene_trees[[a]][[1]][,1] == length(gene_trees[[a]]$tip.label)+1)
ignore_removal <- vector(mode="numeric", length=0)
for (b in 1:length(ignore)){
if (gene_trees[[a]][[1]][,2][[ignore[[b]]]] > length(gene_trees[[a]]$tip.label)){
if (setequal(outgroup_names[[a]], extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,2][[ignore[[b]]]])$tip.label) == FALSE){
ignore_removal <- append(ignore_removal, b)
}
} else {
if (setequal(outgroup_names[[a]], gene_trees[[a]]$tip.label[[gene_trees[[a]][[1]][,2][[ignore[[b]]]]]]) == FALSE){
ignore_removal <- append(ignore_removal, b)
}
}
}
ignore <- ignore[-ignore_removal]
use <- seq(1, nrow(gene_trees[[a]][[1]]), 1)[-ignore]

for (b in use){
if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
gene_trees_descendant <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,2][[b]])
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
} else {
gene_trees_descendant <- gene_trees[[a]]$tip.label[[gene_trees[[a]][[1]][,2][[b]]]]
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
}

if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
species_tree_descendant <- findMRCA(overall_species_tree, gene_trees_descendant$tip.label, "node")
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
} else {
species_tree_descendant <- which(overall_species_tree$tip.label == gene_trees_descendant)
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
}

for (c in 1:nrow(overall_species_tree[[1]])){
if (overall_species_tree[[1]][,2][[c]] == species_tree_descendant & overall_species_tree[[1]][,1][[c]] == species_tree_ancestral){
species_tree_branch_estimate[[c]] <- append(species_tree_branch_estimate[[c]], gene_trees[[a]]$edge.length[[b]])
}
}

}
}

}

########################
###BRANCH_LENGTH_TREE###
########################

branch_length_tree <- overall_species_tree

for (i in 1:length(species_tree_branch_estimate)){
species_tree_branch_estimate[[i]] <- mean(species_tree_branch_estimate[[i]])
}

for (i in 1:length(species_tree_branch_estimate)){
if (is.na(species_tree_branch_estimate[[i]]) == TRUE){
species_tree_branch_estimate[[i]] <- 0
}
}


branch_length_tree$edge.length <- unlist(species_tree_branch_estimate)