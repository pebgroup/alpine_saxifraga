library(ggplot2)
library(sf)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)

coords_and_colours <- fread("longs_lats_colours_two.tsv") 
longitude <- coords_and_colours[[1]]
latitude <- coords_and_colours[[2]]
col <- coords_and_colours[[3]]

points_df <- data.frame(longitude, latitude)
points_df <- st_as_sf(points_df, coords = c("longitude", "latitude"), crs = 4326)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


smaller_points <- st_buffer(points_df, dist = 0.1)

plot_map <- ggplot(data = world) +
geom_sf(data=world, fill = "grey78", colour ="grey78", lty=0)+#, alpha=0)+
geom_sf(data=points_df,size=0.25,stroke=0.25,colour=col,pch=20) +
coord_sf(crs = "+proj=eqearth") + # projection adjusted here
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave("plot_new_colour_friendly3.png",plot_map)

