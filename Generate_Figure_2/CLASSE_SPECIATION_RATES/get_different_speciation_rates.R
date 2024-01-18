library(phytools)
library(ggplot2)

###get_table###

table <- read.table("habitat_parameters.log")

###read in analysis output###

non_alpine <- as.numeric(table[,which(table[1,] == "speciation_in_non_alpine")][-1])
alpine <- as.numeric(table[,which(table[1,] == "speciation_in_alpine")][-1])
both <- as.numeric(table[,which(table[1,] == "speciation_in_both")][-1])

###remove burnin###

burnin = 0.25
non_alpine <- non_alpine[-seq(1, burnin*length(non_alpine), 1)]
alpine <- alpine[-seq(1, burnin*length(alpine), 1)]
both <- both[-seq(1, burnin*length(both), 1)]

###plot_as_histogram###

df <- data.frame(non_alpine, alpine, both)

p <- ggplot(df, aes(x=df[,1])) + 
geom_density(fill="grey", alpha=0.5, weight=0, adjust=2)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line=element_line(size=1,colour="black"), text=element_text(size=16,colour="black"), axis.ticks=element_line(size=1,colour="black"))+
xlim(0, 0.41)+
ylim(0, 40)

