rm (list = ls())
library(tibble)
library(ggplot2)
library(dplyr)
library(stars)
library(ggspatial)
library(showtext)

# The location of all 28 high-quality chicken genomes used in this study ----------
showtext_auto()
cities <- read.csv('import/Denovo_location.CSV',header = T, check.names = F)
cities$Group <- factor(cities$Group,
                       levels = c("1", "2"),
                       labels = c("This study", "Published"))
head(cities)
df_st_as_sf <- st_as_sf(cities,
                        coords = c('longtitude','latitude'), 
                        crs = 4326)
nat.earth <- raster::brick("./NE1_50M_SR_W/NE1_50M_SR_W/NE1_50M_SR_W.tif") 

ggplot(df_st_as_sf)+  
  layer_spatial(nat.earth)+ 
  theme_bw()+
  geom_sf(mapping = aes(fill = Group), shape = 21, size = 2) +
  scale_fill_manual(values = c("#ED2224","#3B54A3")) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  annotation_north_arrow(location = "tl", which_north = F,
                         pad_x = unit(0.03, "in"),  
                         pad_y = unit(0.03, "in"),  
                         style = north_arrow_fancy_orienteering) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        legend.position = 'right',
        legend.title = element_blank(),
        panel.background = element_rect("#f5f6f1"),
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent")
  )+  
  labs(x='', y='',color=NULL);

ggsave("Pan_genome_denovo_samples.pdf", width = 9, height = 4.5)

# The location of 15 Chinese indigenous chicken genomes assembled in this study----------
travel <- read.csv('import/Denovo_location_this_study.CSV',header = T, check.names = F)
travel_st_as_sf <- st_as_sf(travel,
                            coords = c("longtitude", "latitude"),   
                            crs = 4326)
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")[c("adcode", "name", "geometry")]
nat.earth <- raster::brick("./NE1_50M_SR_W/NE1_50M_SR_W/NE1_50M_SR_W.tif")

ggplot(china_map) +
  layer_spatial(nat.earth)+ 
  theme_bw() +
  geom_sf(fill = "white",color='black',size=1.2) +
  geom_sf(data = travel_st_as_sf, shape = 21, fill =  "#ED2224", size = 3) + 
  annotation_north_arrow(location = "tl", which_north = F,    
                         pad_x = unit(0.05, "in"),       
                         pad_y = unit(0.05, "in"),       
                         style = north_arrow_fancy_orienteering)+ 
  theme(axis.text = element_text(size = 12, color = "black"),  
        axis.line = element_blank(),     
        panel.background = element_rect("white"))+
  coord_sf(ylim = c(-3687082,1654989),     
           xlim = c(-3000000,2700000),      
           crs = "+proj=laea +lat_0=40 +lon_0=104")+  
  annotation_scale(location="bl",width_hint=0.3)+ 
  labs(x='', y='',color=NULL)

ggsave("Pan_genome_denovo_this_study.pdf",device = cairo_pdf, width = 7, height = 4.5)

