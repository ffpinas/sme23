
####Francis comment/modified from Gustavs comment####
### These originates from dario and gustav nmi18 work, full script located here M:\marin\nkp16\hob16\arbetskataloger\habitat\Modellering\Script\SGU\1. predictors_prep\sp1577_hob_20180417.R
### help fix all raster issues, and also can create pca images from stack of variables, and finaly can also create a reduced set of PCA images
### final outputs need to be in individual images(geotiff), .
###goal is to have all predictor variables ready for model testing, which is where we decide which ones to include in final model. 
### in the original scrip there are also functions to create a clipped datset of variables - could be good for testing spatial prediction on smaller area...


# install packages ------------------------------------------------------------------------
x <- c("raster", "RStoolbox", "maptools", "caret", "gbm", "rgdal", "doParallel", "tidyverse", "rsample", "terra")
install.packages(x) # warning: uncommenting this may take a number of minutes

#### load packages ####
library(tidyverse)
library(sf)
library(raster)
library(RStoolbox)
library(caret) #wrapper for multiple types of modeling techniques
library(gbm) #package contains BRTs
library(doParallel)
library(readxl)
library(rsample)
library(terra)
#set.seed(10)	

# set paths ------------------------------------------------------------------
getwd()

path2pca = "./inputs/predictors/pca/"
path2mask = "./inputs/mask/"

path_depth_variables = "./inputs/predictors/depth/"
path_Envi_variables = "./inputs/predictors/depth/Envi/"
path_bsEnvi_variables = "./inputs/predictors/bs/"

tevasten_raster_mask = rast(paste(path2mask, "tevasten_raster_mask.tif", sep =""))#using terra
tevasten_raster_mask1 = raster(paste(path2mask, "tevasten_raster_mask.tif", sep =""), CRS("+init=epsg:3006"))#using the raster
crs(tevasten_raster_mask1) <- CRS("+init=epsg:3006")
tevasten_mask_polgyon <- read_sf(paste(path2mask, "tevasten_poly_mask.shp", sep =""))

#raster_depth = rast(paste(path_depth_variables, "by_sme23_Tevasten_1m.tif", sep =""))



#### Functions ---------------------------------------------

# Check predictor rasters

#check function
check = function(x) {
  if(is.na(projection(x))) print("projection NOT set")
  if(projection(x) != projection(mask2)) print (paste(names(x), " different projection", sep=""))
  if(any(dim(x) != dim(mask2))) print (paste(names(x), " different dimensions", sep=""))
  if(extent(x) != extent(mask2)) print(paste(names(x), " different extent", sep=""))}

check2 = function(x) {
  index = integer()
  if(is.na(projection(x))) index[1]=0 else index[1]=1
  if(projection(x) != projection(mask2)) index[2]=0 else index[2]=1
  if(any(dim(x) != dim(mask2))) index[3]=0 else index[3]=1
  if(extent(x) != extent(mask2)) index[4]=0 else index[4]=1
  return(sum(index))
}



### A. compute for various terrain metrics
  # compute for slope, aspect, Topographic Position Index (TPI), Roughness, Terrain Ruggednes Index (TRI)
  # outputs a .tif file
  terra::terrain(tevasten_raster_mask, "slope", filename = paste(path_depth_variables,"/","slope.tif", sep = ""), overwrite = T)
  terra::terrain(tevasten_raster_mask, "aspect", filename = paste(path_depth_variables,"/","aspect.tif", sep = ""), overwrite = T)
  terra::terrain(tevasten_raster_mask, "TPI", filename = paste(path_depth_variables,"/","TPI.tif", sep = ""), overwrite = T) 
  terra::terrain(tevasten_raster_mask, "roughness", filename = paste(path_depth_variables,"/","roughness.tif", sep = ""), overwrite = T)
  terra::terrain(tevasten_raster_mask, "TRI", filename = paste(path_depth_variables,"/","TRI.tif", sep = ""), overwrite = T)
  
### B. create PCA from the obia files of depth raster
    # loop to stack images in all selected groups and turn them into pca, 
    # then reduce to 95 % cum importance and save as tiff
    
  # function to read pca stats (to extract number reduced bands < 95 contribution) in forloop below
    pca_importance <- function(x) {
    vars <- x$sdev^2
    vars <- vars/sum(vars)
    rbind(`Standard deviation` = x$sdev, `Proportion of Variance` = vars, 
          `Cumulative Proportion` = cumsum(vars))
  }

  # read the ENVI extracted rasters from the depth for PCA analysis
  
  Envi_files <- list.files(path_Envi_variables, pattern = "\\.tif$")
  Envi_files
  
  files_stack <- raster::stack()
  for (each_tif_file in Envi_files) {
       raster1 <- raster(paste(path_Envi_variables, each_tif_file, sep = ""))
       raster1 <- projectRaster(raster1, crs = crs(tevasten_raster_mask))
       raster1 <- resample(raster1, tevasten_raster_mask1)
       raster1 <- mask(raster1, tevasten_raster_mask1)
       files_stack = stack(files_stack , raster1)
        }
  names(files_stack)
  
  # run PCA analysis
  pca <- rasterPCA(files_stack, maskCheck=T, spca= T)
  
  saveRSTBX(pca, file = paste(path2pca, "_pca_depth.RData", sep= ""), format = "raster", overwrite=TRUE)
  
  pca$model
  
  # rm(pca)
  #pca <- readRSTBX(paste(path2pca, group, "_pca.RData", sep= ""))
  
  sum_pca <- pca_importance(summary(pca$model))[3,] ## what is the number "3" in this line for?
  
  n_red <- detect_index(sum_pca, function(x) x > 0.95)
  
  x <- raster::stack()
  layer_name <- names(pca$map)
  for (lyr in layer_name){
      newfile <- raster(assign(lyr, pca$map[[lyr]]))
      x <- stack(x, newfile)
    }
  
  pca_red <- dropLayer(x, c((n_red+1):length(sum_pca)))
  
  writeRaster(pca_red, paste(path2pca, "pca_depth.tif", sep= ""), bylayer = T,  overwrite=T, COMPRESS=LZW)

  rm(pca)
  ### C. create PCA from OBIA rasters extracted from "backscatter(bs)" 
  # loop to stack images in all selected groups and turn them into pca, 
  # then reduce to 95 % cum importance and save as tiff

  # read the ENVI extracted rasters from the backscatter for PCA analysis
  
  bsEnvi_files <- list.files(path_bsEnvi_variables, pattern = "\\.tif$")
  bsEnvi_files
  
  files_stack <- raster::stack()
  for (each_tif_file in bsEnvi_files) {
    raster1 <- raster(paste(path_bsEnvi_variables, each_tif_file, sep = ""))
    raster1 <- projectRaster(raster1, crs = crs(tevasten_raster_mask))
    raster1 <- resample(raster1, tevasten_raster_mask1)
    raster1 <- mask(raster1, tevasten_raster_mask1)
    files_stack = stack(files_stack , raster1)
    
  }
  names(files_stack)
  
  
  
  # run PCA
  pca <- rasterPCA(files_stack, maskCheck=T, spca= T)
  
  saveRSTBX(pca, file = paste(path2pca, "_pca_bs.RData", sep= ""), format = "raster", overwrite=TRUE)
  
  pca$model
  
  # rm(pca)
  #pca <- readRSTBX(paste(path2pca, group, "_pca.RData", sep= ""))
  
  sum_pca <- pca_importance(summary(pca$model))[3,] ## what is the number "3" in this line for?
  
  n_red <- detect_index(sum_pca, function(x) x > 0.95)
  
  x <- raster::stack()
  layer_name <- names(pca$map)
  for (lyr in layer_name){
    newfile <- raster(assign(lyr, pca$map[[lyr]]))
    x <- stack(x, newfile)
    rm(newfile)
  }
  
  pca_red <- dropLayer(x, c((n_red+1):length(sum_pca)))
  
  writeRaster(pca_red, paste(path2pca, "pca_bs.tif", sep= ""), bylayer = T,  overwrite=T, COMPRESS=LZW)

  rm(pca)
    
  


#### create XY raster - # From original script -- Another predictor variable? --------------------------------------------------

  
xy <- data.frame(xyFromCell(mask2, 1:ncell(mask2)))

dfx <- xy %>%
  mutate(value = x)

dfx <- rasterFromXYZ(dfx) 

dfx <- dfx * mask2

projection(dfx) <- projection(mask2)

extent(dfx) <- extent(mask2)

dfy <- xy %>%
  mutate(value = y)

dfy <- rasterFromXYZ(dfy) 

dfy <- dfy * mask2

projection(dfy) <- projection(mask2)

extent(dfy) <- extent(mask2)

writeRaster(dfx, paste(way, "xy/", "x.tif", sep= ""),  overwrite=T, COMPRESS=LZW)
writeRaster(dfy, paste(way, "xy/", "y.tif", sep= ""),  overwrite=T, COMPRESS=LZW)


# Extra non related stuff------------------------------------------------------------
resample_mask <- raster("S:/prj/hob16/arbetskataloger/djup_bs_interp/mask_5m.tif")

resample_mask <- disaggregate(resample_mask, 2)

writeRaster(resample_mask, "S:/prj/hob16/arbetskataloger/djup_bs_interp/mask_2_5m.tif", bylayer = T,  overwrite=T, COMPRESS=LZW)


pca_mask <- disaggregate(mask2, 2)

pca_mask <- raster(paste(path2pca, "pca_mask.tif", sep = ""))

writeRaster(pca_mask, paste(path2pca,"pca_mask.tif", sep= ""), bylayer = T,  overwrite=T, COMPRESS=LZW)
