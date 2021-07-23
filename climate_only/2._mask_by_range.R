##############################################################################
### SCRIPT PURPOSE: Clip climate ratser with desired range extent

# Author: Daniel Anstett
# last update:  July 23 2021

## OVERALL WORKFLOW:
# Assumes you have bioclim raster and shape file for the wanted range extent
# Produces new raster that is constrained by range extent shape file
##############################################################################

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

#Import libraries
library(rgdal)
library(sf)
library(tmap)
library(tidyverse)
library(rnaturalearth)
library(raster)
library(cowplot)
library(rgeos)

# Import M.cardinalis range extent as sf polygon
# Used to further constrain size of climate layer
c_range <- st_read("Data/c_range_2.shp")



#Import climate variable, 1981 to 2010 rasters
laea_CMD<-raster("Data/1981_2010/CMD.grd")
laea_MAP<-raster("Data/1981_2010/MAP.grd")
laea_MAT<-raster("Data/1981_2010/MAT.grd")
laea_PAS<-raster("Data/1981_2010/PAS.grd")
laea_EXT<-raster("Data/1981_2010/EXT.grd")

laea_PPT_sm<-raster("Data/1981_2010/PPT_sm.grd")
laea_PPT_wt<-raster("Data/1981_2010/PPT_wt.grd")
laea_Tave_sm<-raster("Data/1981_2010/Tave_sm.grd")
laea_Tave_wt<-raster("Data/1981_2010/Tave_wt.grd")

#re-project all into WGS 1984 (EPSG 4326)
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS
wgs_CMD <- projectRaster(laea_CMD, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_MAP <- projectRaster(laea_MAP, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_MAT <- projectRaster(laea_MAT, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_PAS <- projectRaster(laea_PAS, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_EXT <- projectRaster(laea_EXT, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)

wgs_PPT_sm <- projectRaster(laea_PPT_sm, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_PPT_wt <- projectRaster(laea_PPT_wt, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_Tave_sm <- projectRaster(laea_Tave_sm, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
wgs_Tave_wt <- projectRaster(laea_Tave_wt, crs=EPSG4326) #reproject to WGS 1984 (EPSG 4326)
c_range <- st_transform(c_range, crs = 4326) # reproject to WGS 1984 (EPSG 4326)


#confirm CRS match
crs(c_range) 
crs(wgs_CMD)

# Reduce range extent to exact rectangular extent of M. cardinalis species range
CMD.clip <- raster::crop(wgs_CMD, extent(c_range))
MAP.clip <- raster::crop(wgs_MAP, extent(c_range))
MAT.clip <- raster::crop(wgs_MAT, extent(c_range))
PAS.clip <- raster::crop(wgs_PAS, extent(c_range))
EXT.clip <- raster::crop(wgs_EXT, extent(c_range))

PPT_sm.clip <- raster::crop(wgs_PPT_sm, extent(c_range))
PPT_wt.clip <- raster::crop(wgs_PPT_wt, extent(c_range))
Tave_sm.clip <- raster::crop(wgs_Tave_sm, extent(c_range))
Tave_wt.clip <- raster::crop(wgs_Tave_wt, extent(c_range))

#Write out Raster
writeRaster(CMD.clip, file="Data/CMD.clip.grd", overwrite=TRUE)
writeRaster(MAP.clip, file="Data/MAP.clip.grd", overwrite=TRUE)
writeRaster(MAT.clip, file="Data/MAT.clip.grd", overwrite=TRUE)
writeRaster(PAS.clip, file="Data/PAS.clip.grd", overwrite=TRUE)
writeRaster(EXT.clip, file="Data/EXT.clip.grd", overwrite=TRUE)

writeRaster(PPT_sm.clip, file="Data/PPT_sm.clip.grd", overwrite=TRUE)
writeRaster(PPT_wt.clip, file="Data/PPT_wt.clip.grd", overwrite=TRUE)
writeRaster(Tave_sm.clip, file="Data/Tave_sm.clip.grd", overwrite=TRUE)
writeRaster(Tave_wt.clip, file="Data/Tave_wt.clip.grd", overwrite=TRUE)
