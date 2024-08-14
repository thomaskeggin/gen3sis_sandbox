# This script downloads and subsets the inputs used by Keggin et al. (2023).

# set --------------------------------------------------------------------------
lib <- c("terra","raster","matrixStats","sp","gdistance","geosphere","parallel", "dplyr")
sapply(lib,library, character.only=TRUE)

# load -------------------------------------------------------------------------
load("./keggin_data/temp.Rdata")

depth_rast <-
  rast(geoDepthList)

temp_rast <-
  rast(geoTempList)

# wrangle ----------------------------------------------------------------------
# first 48 time steps
depth_rast <-
  depth_rast[[1:48]]

temp_rast <-
  temp_rast[[1:48]]

# crop to carribbean
car_ext <-
  ext(c(-99.5,-50.5,0.5,39.5))

depth_rast <-
  crop(depth_rast,car_ext)

temp_rast <-
  crop(temp_rast,car_ext)

# get back to making some distance matrices ------------------------------------

# set variables ----------------------------------------------------------------
OutputDir   <- paste0("./input/seascapes/")
crossing_NA <- 0     # Set to 0 (conductance) making land impassible. See gdistance package documentation.
depth_cut   <- -20000    # set the depth cut-off
temp_cut    <- -100    # set lower temperature limit cut-off
# check, or create, output directories -----------------------------------------
# create dirs if they don't exist
if (!dir.exists(paste0(OutputDir,"/distances_full"))){
  dir.create(file.path(OutputDir, "distances_full"))
}

# create landscapes ------------------------------------------------------------
# filter out uninhabitable cells by environmental cut offs
cutTemp  <- temp_rast
cutDepth <- depth_rast

cutDepth[cutDepth < depth_cut] <- NA
cutTemp[cutTemp < temp_cut] <- NA
cutTemp[is.na(cutDepth)] <- NA
cutDepth[is.na(cutTemp)] <- NA

crs(cutDepth) <-
  "+proj=longlat +datum=WGS84"
crs(cutTemp) <-
  "+proj=longlat +datum=WGS84"

# merge all temperature rasters into a single dataframe
masterTemp <-
  as.data.frame(cutTemp, xy=T)

colnames(masterTemp) <- c("x","y",format(round(geoTimes[1:48], 2), nsmall = 2))

# merge all depth rasters into a single dataframe
masterDepth <-
  as.data.frame(cutDepth, xy=T)

colnames(masterDepth) <- c("x","y",format(round(geoTimes[1:48], 2), nsmall = 2))

# explicitly assign rownames
rownames(masterTemp)  <- 1:dim(masterTemp)[1]
rownames(masterDepth) <- 1:dim(masterDepth)[1]

# create and save landscapes object
landscapes <- list(temp = masterTemp, depth = masterDepth)
saveRDS(landscapes, file = file.path(OutputDir, paste0("landscapes.rds",sep="")))

# create distance matrices -----------------------------------------------------
t_start <- dim(depth_rast)[3]  # the starting raster index
t_end   <- 1                     # the final raster index (present day)

for (i in t_start:t_end){
  
  raster_i <- depth_rast[[i]]
  crs(raster_i) <- "+proj=longlat +datum=WGS84"
  age      <- geoTimes[i]
  
  conductObj                     <- raster_i      # this is setting up the conductance (cost of dispersal) values for all marine cells in the raster
  conductObj[!is.na(conductObj)] <- 1             # this gives habitable cells a cost for crossing (1 = no change in cost)
  conductObj[is.na(conductObj)]  <- crossing_NA   # this gives the NA valued cells a cost for crossing (land)
  
  # create a transition object (based on conductance)
  transObj <- transition(raster(conductObj), transitionFunction=min, directions=8) # create matrix with least cost values between each pair of cells (symmetrical?)
  transObj <- geoCorrection(transObj, type = "r", scl = F) * 1000          # correct for map distortion. The output values are in m, the "*1000" converts to km
  # filter by out cells by environmental cut offs
  raster_i        <- mask(raster_i, cutDepth[[i]])                  # filter buy cut offs implemented in the landscapes step
  df_i            <- as.data.frame(raster_i, xy=TRUE, na.rm = TRUE) # this will remove NA cells
  colnames(df_i)  <- c("x","y","depth")
  mat_i_habitable <- data.matrix(df_i)[, 1:2] # convert to matrix of habitable coordinates
  
  # calculate the least-cost distance between points using the transition object and target (habitable) cells
  dist_mat <- costDistance(transObj,
                           mat_i_habitable,
                           mat_i_habitable)
  
  # number rows and columns
  rownames(dist_mat) <- rownames(masterDepth)[which(!is.na(masterDepth[,i+2]))]
  colnames(dist_mat) <- rownames(masterDepth)[which(!is.na(masterDepth[,i+2]))]
  
  # save the distance matrix
  saveRDS(dist_mat,file=file.path(paste0(OutputDir,"/distances_full/distances_full_",i-1,".rds",sep="")))
  
  cat("Done with", round(age, digits = 2), "\n")
  
}
