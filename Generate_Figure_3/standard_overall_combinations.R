library(phangorn)
library(devtools)
library(phytools)


#############
###OPTIONS###
#############

anagenetic_region_gain <- "yes"
anagenetic_biome_gain <- "yes"

anagenetic_region_loss <- "yes" 
anagenetic_biome_loss <- "yes"

extinction_due_to_region_loss <- "yes"
extinction_due_to_biome_loss <- "yes"

biome_vicariance <- "yes"
region_vicariance <- "yes"

biome_dispersal <- "no"
region_dispersal <- "no"

biome_subset_sympatry <- "no"
region_subset_sympatry <- "no"

no_change <- "yes"

########################
###STARTING_VARIABLES###
########################

regions <- 5
biomes <- 2
alpine_region <- "Arctic"

biome_names <- c("Alpine", "Non-alpine")
region_names <- c("Asia", "Europe", "Caucuses", "America", "Arctic")

################################
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

################################
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

#################################
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
