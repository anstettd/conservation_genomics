##############################################################################
### SCRIPT PURPOSE: Find donor sites using only climate

# Author: Daniel Anstett 
# last update:  July 23 2021
##############################################################################

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) 
library(raster) 
library(cowplot)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
#devtools::install_github("ropensci/rnaturalearthhires") # how to install this package
library(rnaturalearthhires)
library(rgeos)
library(RColorBrewer)


### OVERALL WORKFLOW:
## Assumes you have:
# coordinates of focal populations
# future climate projections for populations of conservation,
# climate rasters that have been bounded by range extent of interest
##Produces: 
#Locations that match future climate changes projections for target population
#How those locations match known locations of focal species

##############################################################################

## INPUTS

#Future Climate Data per Restoration Site
gcc.year <- read_csv("Data/timeseries_lat_2 GCMsY_2041_2071.csv")
gcc.season <- read_csv("Data/timeseries_lat_2 GCMsS_2041_2071.csv")
#Selected wanted variables
gcc.year <- gcc.year %>% dplyr::select(GCM:Elevation,
                                       CMD,MAP,MAT,PAS,EXT)
gcc.season <- gcc.season %>% dplyr::select(PPT_sm,PPT_wt,Tave_sm,Tave_wt,)
gcc.clim <- cbind(gcc.year,gcc.season)

#Set up focal 55 populations into sf object
pop_var_raw <- read_csv("Data/paper_ID_site_select.csv")
gen_pop <- pop_var_raw %>% dplyr::select(Long,Lat)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
gen_pop_sf <- st_as_sf(gen_pop,coords=c("Long","Lat"), crs=EPSG4326)
#check data is set up properly
ggplot()+ geom_sf(data = gen_pop_sf)+ ggtitle("Focal Populations")

#Setup sf object for P8
p8_only <- pop_var_raw %>% filter(Paper_ID==8) %>% dplyr::select(Long,Lat)
p8_sf <- st_as_sf(p8_only,coords=c("Long","Lat"), crs=EPSG4326)


# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#st_crs(calo) in WGS 1984


#Climate Raster
# 9 rasters are available but only MAT, MAP and CMD are demonstrated. 
#Annual
MAT.clip <- raster("Data/MAT.clip.grd")
MAP.clip <- raster("Data/MAP.clip.grd")
PAS.clip <- raster("Data/PAS.clip.grd")
EXT.clip <- raster("Data/EXT.clip.grd")
CMD.clip <- raster("Data/CMD.clip.grd")

#Seasonal
PPT_sm.clip <- raster("Data/PPT_sm.clip.grd")
PPT_wt.clip <- raster("Data/PPT_wt.clip.grd")
Tave_sm.clip <- raster("Data/Tave_sm.clip.grd")
Tave_wt.clip <- raster("Data/Tave_wt.clip.grd")


##############################################################################

### Locations that match future climate changes projections for target population
# Worked example with target recipient Site 8

###### Mean Annual Temperature (MAT) ######

# Plot MAT climate Layer with populations and legend
pall_temp<-rev(brewer.pal(n=6, "RdYlBu")) #setup color gradient
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_pop <- tm_shape(MAT.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette=pall_temp)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = F)
clim_MAT_pop
tmap_save(clim_MAT_pop, filename = "Graphs/clim_MAT_pop.pdf",width=5, height=6)

#setup for ssp245 climate change projection for MAT
MAT.S17.ssp245<-MAT.clip #setup object
#Constrain by climate layer by 2040-2070 cliamte change projection +- 1 degree
MAT.S17.ssp245[MAT.clip<gcc.clim$MAT[8]-1] <-NA 
MAT.S17.ssp245[MAT.clip>gcc.clim$MAT[8]+1] <-NA

#setup for ssp555 climate change projection for MAT
MAT.S17.ssp585<-MAT.clip  #setup object
#Constrain by climate layer by 2040-2070 cliamte change projection +- 1 degree
MAT.S17.ssp585[MAT.clip<gcc.clim$MAT[20]-1] <-NA
MAT.S17.ssp585[MAT.clip>gcc.clim$MAT[20]+1] <-NA


#Make with the full climate Layer and focal populations
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_pop <- tm_shape(MAT.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = F)
clim_MAT_pop
tmap_save(clim_MAT_pop, filename = "Graphs/clim_MAT_pop.pdf",width=5, height=6)

#Map only locations matching ssp245 MAT climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_245 <- tm_shape(MAT.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_245
tmap_save(clim_MAT_245, filename = "Graphs/clim_MAT_245.pdf",width=5, height=6)

#Map only locations matching ssp585 MAT climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_MAT_585 <- tm_shape(MAT.S17.ssp585, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAT_585
tmap_save(clim_MAT_585, filename = "Graphs/clim_MAT_585.pdf",width=5, height=6)



###### Mean Annual Precpitation (MAP) ######

# Plot MAP climate Layer with populations and legend
pall_ppt<-brewer.pal(n=5, "YlGnBu") #setup color gradient
tmap_mode("plot")
#tmap_mode("view")
clim_MAP_pop <- tm_shape(MAP.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette=pall_ppt)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = F)
clim_MAP_pop
tmap_save(clim_MAP_pop, filename = "Graphs/clim_MAP_pop.pdf",width=5, height=6)

#setup for ssp245 climate change projection for MAP
MAP.S17.ssp245<-MAP.clip #setup object
#Constrain by climate layer by 2040-2070 cliamte change projection +- 100 mm
MAP.S17.ssp245[MAP.clip<gcc.clim$MAP[8]-100] <-NA
MAP.S17.ssp245[MAP.clip>gcc.clim$MAP[8]+100] <-NA

#setup for ssp585 climate change projection for MAP
MAP.S17.ssp585<-MAP.clip #setup object
#Constrain by climate layer by 2040-2070 cliamte change projection +- 100 mm
MAP.S17.ssp585[MAP.clip<gcc.clim$MAP[20]-100] <-NA
MAP.S17.ssp585[MAP.clip>gcc.clim$MAP[20]+100] <-NA


#Make with the full climate Layer and focal populations
tmap_mode("plot")
#tmap_mode("view")
clim_MAP_pop <- tm_shape(MAP.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = F)
clim_MAP_pop
tmap_save(clim_MAP_pop, filename = "Graphs/clim_MAP_pop.pdf",width=5, height=6)

#Map only locations matching ssp245 MAP climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_MAP_245 <- tm_shape(MAP.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAP_245
tmap_save(clim_MAP_245, filename = "Graphs/clim_MAP_245.pdf",width=5, height=6)

#Map only locations matching ssp585 MAT climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_MAP_585 <- tm_shape(MAP.S17.ssp585, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_MAP_585
tmap_save(clim_MAP_585, filename = "Graphs/clim_MAP_585.pdf",width=5, height=6)



###### Cumulative Moisture Deficite (CMD) ######

# Plot CMD climate layer with populations and legend
tmap_mode("plot")
#tmap_mode("view")
clim_CMD_pop <- tm_shape(CMD.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = T)
clim_CMD_pop
tmap_save(clim_CMD_pop, filename = "Graphs/clim_CMD_pop.pdf",width=5, height=6)

#setup for ssp245 climate change projection for CMD
CMD.S17.ssp245<-CMD.clip
#Constrain by climate layer by 2040-2070 climate change projection +- 100 units
CMD.S17.ssp245[CMD.clip<gcc.clim$CMD[8]-100] <-NA
CMD.S17.ssp245[CMD.clip>gcc.clim$CMD[8]+100] <-NA

#Constrain by climate layer by 2040-2070 climate change projection +- 100 units
CMD.S17.ssp585<-CMD.clip
CMD.S17.ssp585[CMD.clip<gcc.clim$CMD[20]-100] <-NA
CMD.S17.ssp585[CMD.clip>gcc.clim$CMD[20]+100] <-NA


#Make with the full climate Layer and focal populations
tmap_mode("plot")
#tmap_mode("view")
clim_CMD_pop <- tm_shape(CMD.clip, bbox=st_bbox(calo)) + #legal boundires
  tm_raster()+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_layout(legend.show = F)
clim_CMD_pop
tmap_save(clim_CMD_pop, filename = "Graphs/clim_CMD_pop.pdf",width=5, height=6)

#Map only locations matching ssp245 CMD climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_CMD_245 <- tm_shape(CMD.S17.ssp245, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_CMD_245
tmap_save(clim_CMD_245, filename = "Graphs/clim_CMD_245.pdf",width=5, height=6)

#Map only locations matching ssp585 CMD climate change projections
tmap_mode("plot")
#tmap_mode("view")
clim_CMD_585 <- tm_shape(CMD.S17.ssp585, bbox=st_bbox(calo)) + #legal boundires
  tm_raster(palette = "#FF00FF",legend.show = FALSE)+
  tm_shape(calo)+
  tm_borders()+
  tm_shape(gen_pop_sf)+
  tm_dots(size=0.2,shape=1)+
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
clim_CMD_585
tmap_save(clim_CMD_585, filename = "Graphs/clim_CMD_585.pdf",width=5, height=6)


