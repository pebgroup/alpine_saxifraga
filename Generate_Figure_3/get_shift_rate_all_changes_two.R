library(phytools)
library(NELSI)
library(ggplot2)

###initial_info###

###tree

index_tree <- read.annotated.nexus("index.tre")

branch_indices <- vector(mode="numeric", length=0)
for (i in 1:length(index_tree[[6]])){
branch_indices <- append(branch_indices, index_tree[[6]][[i]]$index)
}

n_branches <- length(index_tree$tip.label) * 2 - 2
table <- data.frame(index_tree[[1]], seq(1, n_branches, 1), branch_indices[-length(branch_indices)], index_tree$edge.length)
table <- table[order(table[,2]),] 
end_times <- c(rep(0, length(index_tree$tip.label)), branching.times(index_tree)[-1])
start_times <- end_times + table[,5]
table <- cbind(table, end_times, start_times)

###biomes

source("standard_overall_combinations.R")
alpine_states <- vector(mode="numeric", length=0)
all_alpine <- vector(mode="numeric", length=0)
exclusively_alpine <- vector(mode="numeric", length=0)
for (a in 1:length(overall_combinations_list)){
alpine_counter <- vector(mode="numeric", length=0)
biome_counter <- vector(mode="numeric", length=0)
for (b in 1:length(overall_combinations_list[[a]])){
if (overall_combinations_list[[a]][[b]] == "Alpine"){
alpine_counter <- append(alpine_counter, 1)
}
if (overall_combinations_list[[a]][[b]] %in% biome_names){
biome_counter <- append(biome_counter, 1)
}
}
if ((length(biome_counter) == 1) & (length(alpine_counter) == 1)){
exclusively_alpine <- append(exclusively_alpine, a)
}
if (((length(biome_counter) == 1) & (length(alpine_counter) == 1)) || (length(biome_counter) == 2)){
all_alpine <- append(all_alpine, a)
}
}

##################

info <- read.table("habitat_smap_clean.log")
info <- info[,-1]

rep_changes_accross <- vector("list", ((nrow(info)-round(0.25*nrow(info)))+1)) # i.e a list the length of the number of generations
rep_times_accross <- vector("list", ((nrow(info)-round(0.25*nrow(info)))+1))

rep_changes_into <- vector("list", ((nrow(info)-round(0.25*nrow(info)))+1))
rep_times_into <- vector("list", ((nrow(info)-round(0.25*nrow(info)))+1))

for (a in round(0.25 * nrow(info)):nrow(info)){
print(a)
changes_into <- vector("list", 0) # revectoring for each generation
times_into <- vector(mode="numeric", length=0)
changes_accross <- vector("list", 0)
times_accross <- vector(mode="numeric", length=0)
for (b in 1:(ncol(info)-2)){ # -2 is because we are not interested in the ancestor to the root node
if (length(unlist(strsplit(as.character(info[a,][b]), split=",")))>2){ #if there are more than two commas there is a change
sections <- unlist(strsplit(as.character(info[a,][b]), split=":")) # extracting the section corresponding to each state
counter <- vector(mode="numeric", length=0)
for (c in 1:(length(sections)-1)){
counter <- append(counter, as.numeric(unlist(strsplit(sections[[c]], split=","))[[2]])) # cumulative ages from which to calculate the correct time 
if ((as.numeric(unlist(strsplit(sections[[c+1]], split=","))[[1]]) %in% exclusively_alpine == FALSE) & (as.numeric(unlist(strsplit(sections[[c]], split=","))[[1]]) %in% exclusively_alpine == TRUE)){
changes_into[[length(changes_into)+1]] <- c(unlist(strsplit(sections[[c+1]], split=","))[[1]], unlist(strsplit(sections[[c]], split=","))[[1]])
times_into <- append(times_into, as.numeric(table[,6][which(table[,4] == b)]) + sum(counter)) # b starts at 1
}
if ((as.numeric(unlist(strsplit(sections[[c+1]], split=","))[[1]]) %in% exclusively_alpine == TRUE) & (as.numeric(unlist(strsplit(sections[[c]], split=","))[[1]]) %in% exclusively_alpine == TRUE)){
if (as.numeric(unlist(strsplit(sections[[c+1]], split=","))[[1]]) < as.numeric(unlist(strsplit(sections[[c]], split=","))[[1]])){ 
changes_accross[[length(changes_accross)+1]] <- c(unlist(strsplit(sections[[c+1]], split=","))[[1]], unlist(strsplit(sections[[c]], split=","))[[1]])
times_accross <- append(times_accross, as.numeric(table[,6][which(table[,4] == b)]) + sum(counter)) 
}
}
}
final_state <- unlist(strsplit(sections[[length(sections)]], split=","))[[1]]
if (table[,1][which(table[,4] == b)] != length(index_tree$tip.label) + 1){##need filter on it not being a root branch
ancestral_state <- unlist(strsplit(as.character(info[a,][table[,4][which(table[,2] == table[,1][which(table[,4] == b)])]]), split=","))[[1]]
} else {
ancestral_state <- unlist(strsplit(as.character(info[a,][length(info[a,])-1]), split=","))[[1]]
}
if (final_state != ancestral_state){
if ((as.numeric(ancestral_state) %in% exclusively_alpine == FALSE) & (as.numeric(final_state) %in% exclusively_alpine == TRUE)){
changes_into[[length(changes_into)+1]] <- c(ancestral_state, final_state)
times_into <- append(times_into, as.numeric(table[,6][which(table[,2] == table[,1][which(table[,4] == b)])]))
}
if ((as.numeric(ancestral_state) %in% exclusively_alpine == TRUE) & (as.numeric(final_state) %in% exclusively_alpine == TRUE)){
if (as.numeric(ancestral_state) < as.numeric(final_state)){
changes_accross[[length(changes_accross)+1]] <- c(ancestral_state, final_state)
times_accross <- append(times_accross, as.numeric(table[,6][which(table[,2] == table[,1][which(table[,4] == b)])]))
}
}
}
} else {
final_state <- unlist(strsplit(as.character(info[a,][b]), split=","))[[1]]
if (table[,1][which(table[,4] == b)] != length(index_tree$tip.label) + 1){
ancestral_state <- unlist(strsplit(as.character(info[a,][table[,4][which(table[,2] == table[,1][which(table[,4] == b)])]]), split=","))[[1]] 
} else {
ancestral_state <- unlist(strsplit(as.character(info[a,][length(info[a,])-1]), split=","))[[1]] 
}
if (final_state != ancestral_state){
if ((as.numeric(ancestral_state) %in% exclusively_alpine == FALSE) & (as.numeric(final_state) %in% exclusively_alpine == TRUE)){
changes_into[[length(changes_into)+1]] <- c(ancestral_state, final_state)
times_into <- append(times_into, as.numeric(table[,6][which(table[,2] == table[,1][which(table[,4] == b)])]))
}
if ((as.numeric(ancestral_state) %in% exclusively_alpine == TRUE) & (as.numeric(final_state) %in% exclusively_alpine == TRUE)){
if (as.numeric(ancestral_state) < as.numeric(final_state)){
changes_accross[[length(changes_accross)+1]] <- c(ancestral_state, final_state)
times_accross <- append(times_accross, as.numeric(table[,6][which(table[,2] == table[,1][which(table[,4] == b)])]))
}
}
}
}
rep_changes_accross[[a]] <- changes_accross ###need to add extension
rep_times_accross[[a]] <- times_accross
rep_changes_into[[a]] <- changes_into
rep_times_into[[a]] <- times_into
}
}

rep_changes_accross <- rep_changes_accross[-seq(1, 138, 1)]
rep_times_accross <- rep_times_accross[-seq(1, 138, 1)]
rep_changes_into <- rep_changes_into[-seq(1, 138, 1)]
rep_times_into <- rep_times_into[-seq(1, 138, 1)]

##put_in_time_windows###
########################

old_into <- vector(mode="numeric", length=0)
old_accross <- vector(mode="numeric", length=0)

into_15_10 <- vector(mode="numeric", length=0)
accross_15_10 <- vector(mode="numeric", length=0)

into_10_5 <- vector(mode="numeric", length=0)
accross_10_5 <- vector(mode="numeric", length=0)

young_into <- vector(mode="numeric", length=0)
young_accross <- vector(mode="numeric", length=0)

all_into <- vector(mode="numeric", length=0)
all_accross <- vector(mode="numeric", length=0)

###

for (a in 1:length(rep_changes_accross)){

old_into <- append(old_into, length(which(rep_times_into[[a]] > 15)))
old_accross <- append(old_accross, length(which(rep_times_accross[[a]] > 15)))

into_15_10 <- append(into_15_10, length(which(rep_times_into[[a]] <= 15 & rep_times_into[[a]] > 10)))
accross_15_10 <- append(accross_15_10, length(which(rep_times_accross[[a]] <= 15 & rep_times_accross[[a]] > 10)))

into_10_5 <- append(into_10_5, length(which(rep_times_into[[a]] <= 10 & rep_times_into[[a]] > 5)))
accross_10_5 <- append(accross_10_5, length(which(rep_times_accross[[a]] <= 10 & rep_times_accross[[a]] > 5)))

young_into <- append(young_into, length(which(rep_times_into[[a]] <= 5)))
young_accross <- append(young_accross, length(which(rep_times_accross[[a]] <= 5)))

all_into <- append(all_into, length(rep_times_into[[a]]))
all_accross <- append(all_accross, length(rep_times_accross[[a]]))

}

##get_branch_length_in_time_window##
####################################

end = 6
start =7

lengths_all <- sum(index_tree$edge.length)
lengths_15_10 <- vector(mode="numeric", length=0)
lengths_10_5 <- vector(mode="numeric", length=0)
lengths_5_0 <- vector(mode="numeric", length=0)

###

for (a in 1:nrow(table)){
if (table[,7][[a]] <= 15 & table[,6][[a]] > 10){
lengths_15_10 <- append(lengths_15_10, table[,5][[a]])
} else if (table[,7][[a]] > 15 & table[,6][[a]] < 10){  ## overshoots either end
lengths_15_10 <- append(lengths_15_10, 5)
} else if (table[,7][[a]] > 15 & table[,6][[a]] <= 15 & table[,6][[a]] > 10){  ## overshoots top
lengths_15_10 <- append(lengths_15_10, 15-table[,6][[a]])
} else if (table[,7][[a]] <= 15 & table[,7][[a]] > 10 & table[,6][[a]] < 10){ ## overshoots bottom
lengths_15_10 <- append(lengths_15_10, table[,7][[a]]-10)
}
}

###

for (a in 1:nrow(table)){
if (table[,7][[a]] <= 10 & table[,6][[a]] > 5){
lengths_10_5 <- append(lengths_10_5, table[,5][[a]])
} else if (table[,7][[a]] > 10 & table[,6][[a]] < 5){  ## overshoots either end
lengths_10_5 <- append(lengths_10_5, 5)
} else if (table[,7][[a]] > 10 & table[,6][[a]] <= 10 & table[,6][[a]] > 5){  ## overshoots top
lengths_10_5 <- append(lengths_10_5, 10-table[,6][[a]])
} else if (table[,7][[a]] <= 10 & table[,7][[a]] > 5 & table[,6][[a]] < 5){ ## overshoots bottom
lengths_10_5 <- append(lengths_10_5, table[,7][[a]]-5)
}
}

###

for (a in 1:nrow(table)){
if (table[,7][[a]] <= 5 & table[,6][[a]] >= 0){
lengths_5_0 <- append(lengths_5_0, table[,5][[a]])
} else if (table[,7][[a]] > 5 & table[,6][[a]] <= 5 & table[,6][[a]] >= 0){  ## overshoots top
lengths_5_0 <- append(lengths_5_0, 5-table[,6][[a]])
}
}

###

lengths_15_and_above <- lengths_all - sum(lengths_15_10) - sum(lengths_10_5) - sum(lengths_5_0)

###

lengths_all <- sum(index_tree$edge.length)

###get_vectors###

density_old_into <- old_into/lengths_15_and_above
density_old_accross <- old_accross/lengths_15_and_above

density_15_10_into <- into_15_10/sum(lengths_15_10)
density_15_10_accross <- accross_15_10/sum(lengths_15_10)

density_10_5_into <- into_10_5/sum(lengths_10_5)
density_10_5_accross <- accross_10_5/sum(lengths_10_5)

density_5_0_into <- young_into/sum(lengths_5_0)
density_5_0_accross <- young_accross/sum(lengths_5_0)

density_into <- all_into/lengths_all
density_accross <- all_accross/lengths_all

###plot

starts <- c(80, 15, 10, 5)
ends <- c(15.0001, 10.00001, 5.00001, 0)

time_intervals <- vector(mode="numeric", length=0)
for (i in 1:length(starts)){
time_intervals <- append(time_intervals, starts[[i]])
time_intervals <- append(time_intervals, ends[[i]]+0.0001)
}

into_means <- rep(c(mean(density_old_into), mean(density_15_10_into), mean(density_10_5_into), mean(density_5_0_into)), each=2)
accross_means <- rep(c(mean(density_old_accross), mean(density_15_10_accross), mean(density_10_5_accross), mean(density_5_0_accross)), each=2)

###through time###

rolling_table <- data.frame(time_intervals, into_means, accross_means) 

rolling_plot <- ggplot(data=rolling_table, aes(x = rolling_table[,1], y = rolling_table[,2]))+
geom_line(aes(y=rolling_table[,2]), colour="black", size=1.2) + 
geom_line(aes(y=rolling_table[,3]), colour="blue", size=1.2) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=16,colour="black"), axis.ticks=element_line(size=1,colour="black"))+
labs(y="rate", x="Ma")+
ylim(0, 0.08)+
scale_x_reverse(limits=c(80,0))

###density###

df <- data.frame(density_accross)

plot <- ggplot(df, aes(x=df[,1])) + 
geom_density(fill="grey", alpha=0.5, weight=0, adjust=2)+
#geom_segment(x=overall_diversification_in_alpine_lower, xend=overall_diversification_in_alpine_lower, y=11, yend=12)+ 
#geom_segment(x=overall_diversification_in_alpine_upper, xend=overall_diversification_in_alpine_upper, y=11, yend=12)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=16,colour="black"), axis.ticks=element_line(size=1,colour="black"))+
xlim(0, 0.125)+
ylim(0, 75)



















