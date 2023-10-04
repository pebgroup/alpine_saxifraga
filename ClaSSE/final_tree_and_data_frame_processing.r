library(phangorn)
library(devtools)
library(phytools)
library(phylobase)

###########################################
###GENERATE_CHARACTER_STATE_COMBINATIONS###
###########################################

###STARTING_VARIABLES###
########################

regions <- 5
biomes <- 2
alpine_region <- "Arctic"

biome_names <- c("Alpine", "Non-alpine")
region_names <- c("Asia", "Europe", "Caucuses", "America", "Arctic")

###GET_COMBINATIONS_OF_REIONS###
################################

region_combinations <- vector("list", 0) 
for (i in 1:regions){
region_combinations[[i]] <- combn(region_names, i)
}

region_combinations_list <- vector("list", 0)
for (j in 1:length(region_combinations)){
for (i in 1:ncol(region_combinations[[j]])){
region_combinations_list[[length(region_combinations_list)+1]] <- region_combinations[[j]][,i]
}
}

###GET_COMBINATIONS_OF_BIOMES###
################################

biome_combinations <- vector("list", 0) 
for (i in 1:biomes){
biome_combinations[[i]] <- combn(biome_names, i)
}

biome_combinations_list <- vector("list", 0)
for (j in 1:length(biome_combinations)){
for (i in 1:ncol(biome_combinations[[j]])){
biome_combinations_list[[length(biome_combinations_list)+1]] <- biome_combinations[[j]][,i]
}
}

###GET_OVERALL_COMBINATION_SET###
#################################

overall_combinations_list <- rep(region_combinations_list, each=length(biome_combinations_list))
biome_combinations_list_storage <- biome_combinations_list
biome_combinations_list <- rep(biome_combinations_list, length(overall_combinations_list)/length(biome_combinations_list))
for (i in 1:length(overall_combinations_list)){
overall_combinations_list[[i]] <- c(overall_combinations_list[[i]], biome_combinations_list[[i]])
}

removal <- vector(mode="numeric", length=0)
for (i in 1:length(overall_combinations_list)){
if (length(which(region_names %in% overall_combinations_list[[i]])) == 1){
if (region_names[[which(region_names %in% overall_combinations_list[[i]])]] == alpine_region){
if ("Non-alpine" %in% overall_combinations_list[[i]]){
removal <- append(removal, i)
}
}
}
if (length(which(biome_names %in% overall_combinations_list[[i]])) == 1){
if (biome_names[[which(biome_names %in% overall_combinations_list[[i]])]] == "Non-alpine"){
if (alpine_region %in% overall_combinations_list[[i]]){
removal <- append(removal, i)
}
}
}
}
removal <- unique(removal)

overall_combinations_list <- overall_combinations_list[-removal]

##########################################
###MATCH_TIP_MATRIX_TO_CHARACTER_STATES###
##########################################

###MATCHING_UP###
#################

matrix_order <- c("Asia", "Europe", "Caucuses", "America", "Arctic", "Alpine", "Non-alpine")
individualised_matrix <- read.csv("biome_region.csv", header=FALSE)

final_matrix_value <- vector(mode="numeric", length=0)
for (j in 1:nrow(individualised_matrix)){
temp_record <- vector(mode="numeric", length=0)
for (i in 2:length(individualised_matrix[j,])){
if (individualised_matrix[j,][[i]] == 1){
temp_record <- append(temp_record, matrix_order[[i-1]])
}
}
for (i in 1:length(overall_combinations_list)){
if (setequal(sort(overall_combinations_list[[i]]), sort(temp_record)) == TRUE){
final_matrix_value <- append(final_matrix_value, i)
}
}
}

###CHARACTER_STATE_TABLE###
###########################

final_table <- data.frame(individualised_matrix[,1], final_matrix_value)

###########################
###PROCESS_TREE_TO_TABLE###
###########################

###READ_IN_TREE###
##################

tree <- read.tree("...") ### need output tree from treePL to prune for subsequent classe analysis
tree <- drop.tip(tree, "Saxifraga_crustata_P07_WC05")

for (i in 1:length(tree$tip.label)){
tree$tip.label[[i]] <- paste(unlist(strsplit(tree$tip.label[[i]][[1]], split="_"))[[1]], unlist(strsplit(tree$tip.label[[i]][[1]], split="_"))[[2]], sep="_")
}

###EXTRACT_SAXAFRAGA_CLADE_AND_REMOVE_DUP_TIPS###
#################################################

tree <- extract.clade(tree, findMRCA(tree, c("Cascadia_nuttallii", "Saxifraga_thianthia"), "node"))

tip_storage <- vector(mode="numeric", length=0)
tip_removal_vector <- vector(mode="numeric", length=0)
for (i in 1:length(tree$tip.label)){
if ((tree$tip.label[[i]] %in% tip_storage) == FALSE){ 
tip_storage <- append(tip_storage, tree$tip.label[[i]]) 
} else if ((tree$tip.label[[i]] %in% tip_storage) == TRUE){ 
tip_removal_vector <- append(tip_removal_vector, i)
}
}

one_per_taxon_tree <- drop.tip(tree, tip_removal_vector)

####################
###SORT_OUT_TABLE###
####################

#after hyphens removed from table
#hyphens amalgamated for matta
#styriaca corrected in tree
#poluniana corrected in table 

missing_table_names <- vector(mode="numeric", length=0)
for (i in 1:nrow(final_table)){
if ((final_table[,1][[i]] %in% one_per_taxon_tree$tip.label) == FALSE){
missing_table_names <- append(missing_table_names,  i)
}
}

###################
###FINALISE_TREE###
###################
final_drop_vector <- vector(mode="numeric", length=0)
if (length(missing_table_names) == 0){
for (i in 1:length(one_per_taxon_tree$tip.label)){
if ((one_per_taxon_tree$tip.label[[i]] %in% final_table[,1]) == FALSE){
final_drop_vector <- append(final_drop_vector, i)
}
}
print("ABLE_TO_RUN")
}

if (length(final_drop_vector)>0){
final_tree <- drop.tip(one_per_taxon_tree, final_drop_vector)
} else {
final_tree <- one_per_taxon_tree
}

##################
###WRITE_OUTPUT###
##################

write.tree(final_tree, "analysis_tree.tre")
write.table(final_table, "analysis_matrix.csv", sep=",", col.names=FALSE, row.names=FALSE)





























