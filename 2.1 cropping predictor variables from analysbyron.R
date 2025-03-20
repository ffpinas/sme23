x <- c("raster", "RStoolbox", "maptools", "caret", "gbm", "rgdal", "doParallel", "tidyverse", "rsample", "terra")
install.packages(x) # warning: uncommenting this may take a number of minutes

#### load packages ####
library(tidyverse)
library(sf)
library(raster)
library(RStoolbox)
library(caret) #wrapper for multiple types of modeling techniques
library(gbm) #package contains BRTs
#library(rgdal) #package lets you I/O rasters
library(doParallel)
library(readxl)
library(rsample)
library(terra)

#script to crop Oscar's predictor variables using prepared mask--------------------------------------------
getwd()


path2pred = "./inputs/predictors/"
path2mask = "./inputs/mask/"


raster_mask = "inputs/mask/tevasten_raster_mask.tif"
polygon_mask <- "inputs/mask/tevasten_poly_mask.shp"


tevasten_raster_mask = rast(paste(path2mask, "tevasten_raster_mask.tif", sep =""))#using terra
tevasten_raster_mask1 = raster(paste(path2mask, "tevasten_raster_mask.tif", sep =""))#using the raster
tevasten_polygon_mask <- read_sf(paste(path2mask, "lovelse_poly_mask.shp", sep =""))



#st_layers(mask_path) # to check the name of the mask file
#tevasten_mask <- read_sf(mask_path, layer = "lovelse_poly_mask" #import the polygon
#extent(lovelse_mask)

#raster mask (use the depth)
tevasten_mask_depth <- rast(raster_mask)
tevasten_mask_polgyon <- read_sf(polygon_mask)

x <- crs(tevasten_raster_mask)

#list all the files to crop
path2pred_Oscar <- "M:/marin/nkp16/sme23/arbetskataloger/Oscar_Baltic100M_analybyran/Batic_100m"

file_list <- grep("\\d", list.files(path2pred_Oscar, pattern = ".tif$", recursive = T, full.names = T), value = T)
Oscar_predictors <- grep("old", file_list, invert = T, value = T)#remove tiff files with "old"
Oscar_predictors1 <- c(Oscar_predictors[1],Oscar_predictors[2], Oscar_predictors[3])

for(raster_files in Oscar_predictors)  {
    folder_name <- basename(dirname(raster_files))
    name <- str_sub((basename(raster_files)), end = -5)
    raster1 <- raster(raster_files)
    raster1 <-projectRaster(raster1, crs = crs(tevasten_raster_mask))
    raster1 <- terra::rast(raster1)
    raster1 <- terra::focal(raster1, w=5, fun=mean, NAonly=T, na.rm=T)
    raster1 <- terra::resample(raster1, tevasten_raster_mask)
    raster1 <- terra::crop(raster1, tevasten_raster_mask, mask=T)
    writeRaster(raster1, paste(path2pred, folder_name,"/", name, "1.tif", sep = ""), overwrite=T)
    
    }


#Use focalvalues to expand the rasters

# Apply the focal function with a sum operation

result <- focalValues(raster_cp, w=3, row=1, nrows=nrow(x), fill=NA)
raster_cp <- rast(raster_cp)
