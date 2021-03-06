##############################################################################
### SCRIPT PURPOSE: Transform climate rasters to bioclim 
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert & Daniel Anstett
# last update:  July 23 2021

## OVERALL WORKFLOW:
# !!! For raster data, this assumes you have downloaded tiff rasters for 36 annual variables 
#from climateNA (https://adaptwest.databasin.org/pages/adaptwest-climatena)
##############################################################################


############################################################################## 
### LOAD LIBRARIES AND PREPARE INPUTS

## CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))

## LIBRARIES
library(tidyverse) # for data manipulation
library(dismo) # for biovars function
library(raster) # for raster grids
library(rgdal) # for transforming projections

## INPUTS
## points defining range of M. cardinalis for the sake of cropping climate layers
 clim <- read_csv("Data/points_Normal_1961_1990MSY.csv")

# Ensure no missing values
complete <- clim[complete.cases(clim), ] # should be same dim as clim if there are no missing values

##############################################################################
#Functions

# A function to crop climate rasters by a bounding box
clip <- function(raster, shape) { 
  a1_crop <- crop(raster, shape)
  step1 <- rasterize(shape, a1_crop)
  return(a1_crop*step1)
}

##############################################################################
### MANIPULATE RASTER FILES TO DECREASE SIZE TO RELEVANT AREA

## Read in TIF files all at once
# !!! Set working directory to folder with ASCII files downloaded from https://www.dropbox.com/sh/thd90znbkbmdfio/AACJfYY_7QLBvvUikiEfZKjWa?dl=0
setwd("~/Dropbox/AM_Workshop/Climate/Normal_1981_2010_bioclim/Selected")
allfiles.list <- list.files(pattern = '.tif') # list of '.tif' files
allfiles <- stack(allfiles.list) # import set of rasters
plot(allfiles[[2]]) # optional visual check, plotting "plot second raster"

#CRS manipulation
crs(allfiles[[1]]) # get CRS. In laea projection, Lambert azimuthal equal-area projection
prj.laea <- "+proj=laea +lon_0=-100 +lat_0=45 +type=crs" ## set laea CRS
proj4string(allfiles) <- CRS(prj.laea) #define projection of rasters
projection(allfiles) #optional check correct projection applied

# The rasters are for all of North America and too large for storage in this repo
# Trim them to study area and save only trimmed files in the repo
# Define extent as 1 degree beyond lat-long extent of points
ext <- extent(min(clim$Longitude)-1, max(clim$Longitude)+1, min(clim$Latitude)-1, max(clim$Latitude)+1)
bbox = as(ext, "SpatialPolygons") #convert coordinates to a bounding box
prj.wgs = "+proj=longlat + type=crs" #define unprojected coordinate system
proj4string(bbox) <- CRS(prj.wgs) #set projection
bbox.lcc = spTransform(bbox, CRS=CRS(prj.laea)) #re-project to match rasters

#mask(i.e. clip) files by bounding box
allfiles.clip <- clip(allfiles, bbox.lcc) #apply the function to the raster stack
names(allfiles.clip) <- names(allfiles) #replace layer names
plot(allfiles.clip[[2]]) #plot clipped layer (PPT02 as an example)
plot(allfiles[[2]]) #compare to unclipped layer

##############################################################################
# Save individual rasters
# !!! Be sure working directory is set back to project folder
setwd("~/Dropbox/a Papers/final_githubs/conservation_genomics")

writeRaster(allfiles.clip[[1]],"Data/1981_2010/CMD.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[2]],"Data/1981_2010/EXT.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[3]],"Data/1981_2010/MAP.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[4]],"Data/1981_2010/MAT.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[5]],"Data/1981_2010/PAS.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)

writeRaster(allfiles.clip[[6]],"Data/1981_2010/PPT_sm.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[7]],"Data/1981_2010/PPT_wt.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[8]],"Data/1981_2010/Tave_sm.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
writeRaster(allfiles.clip[[9]],"Data/1981_2010/Tave_wt.grd",bylayer = TRUE,
            datatype = 'INT4S', bandorder = '.grd', overwrite = TRUE)
##############################################################################
