library(rgdal)
library (raster)
library(sf)
library(dplyr)

sat = brick("satellite_data/Lsat_q1_2022_Jan_GLCMofBands_and_indices_terra.tif")

# Import a polygon shapefile
leakage = st_read("data2/shapefiles/Leakage_Belt_clipped_to_1kmbufferEnvelope_for_Predictions.shp")
grnp = st_read("data2/shapefiles/GRNP_clipped_to_1kmBuffer_Envelope_for_Predictions.shp")

# leakage = st_transform(leakage, CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))
# grnp = st_transform(leakage, CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))
plot(sat[[1]])
plot(grnp$geometry, add=T)
plot(leakage$geometry, add=T)

######
## To Get Distance From GRNP --- I Need to create 2 data sets
## Points INSIDE & Points OUTSIDE
## Then when calculating Distance I make the inside points -ve values

# This creates a grid & points over an sf object (Such as a landsat image, or the Bounding Box around it)
grid_BUFF = st_make_grid(x = sat[[1]], cellsize = 150, what = "centers") %>%
  st_intersection(leakage$geometry)
coords = st_coordinates(grid_BUFF)
coords = as.data.frame(coords[,c(1:2)])
coords$grid_number = rownames(coords)
# write.csv(coords, file = "../data2/prediction_coordinate_150mgrid_gola_buffer_EX_GIS.csv", row.names = F)
xy.out = coords[,1:2]  

grid_GRNP = st_make_grid(x = sat[[1]], cellsize = 150, what = "centers") %>%
  st_intersection(grnp$geometry)
coords = st_coordinates(grid_GRNP)
coords = as.data.frame(coords[,c(1:2)])
coords$grid_number = rownames(coords)
coords$grid_number = paste0("inside_", coords$grid_number)
# write.csv(coords, file = "data2/prediction_coordinate_150mgrid_gola_buffer_EX_GIS.csv", row.names = F)
xy.in = coords[,1:2]

### NB: BOTH SETS OF COORDINATES SUBSEQUENTLY GET DISTANCE TO PARK EDGE & ARE CONVERTED TO PROJECTED CRS IN ArcGIS Pro

#######################################################################################

# For each set of coordinates get buffer data
# run below for 1 then the other
xy = xy.in
# THEN RUN
# xy = xy.out

spdf <- SpatialPointsDataFrame(coords = xy, data = coords,
                               proj4string = CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))

# # ----------------------------
# # Import the satellite images
# # ----------------------------
rad = brick("data2/Sent1_2022_q1_Radar.tif")
hansen = raster("data2/Hansen_GFC-2019-v1.7_lossyear_10N_020W.tif")

### CREATE VALUES FOR EACH BUFFER RADIUS ###
# 1. Radar - Sentinel-1 
p = spTransform(spdf, crs(rad))
# plot(rad[[1]])
# plot(p, add=T)

conflict_prefer("extract", "raster")

rad_mean = extract(x = rad,
                   y = p,
                   buffer=30,
                   fun=mean,
                   df=F)
rad_SD = extract(x = rad,
                   y = p,
                   buffer=30,
                   fun=sd,
                   df=F)
radar = cbind(rad_mean, rad_SD)
colnames(radar) = c("VV_mean", "VH_mean", "VV_sd", "VH_sd")
colnames(radar) = paste0(colnames(radar), "_30")
# file = paste0("data2/Prediction_150mGrid_GRNP_Sentinel1_2022_q1_30m.csv")
# file = paste0("data2/Prediction_150mGrid_buffer_Sentinel1_2022_q1_30m.csv")
# write.csv(x = radar, file = file, row.names = F)

# 2. Hansen year of deforesetation - NB I use Max & Min of year within the buffer
p = spTransform(spdf, crs(hansen))
hansYr_max = extract(x = hansen,
                     y = p,
                     buffer = 30,
                     fun = max,
                     df=F)
  
hansYr_min = extract(x = hansen,
                     y = p,
                     buffer = 30,
                     fun = min,
                     df=F)
hansYr = cbind(hansYr_max, hansYr_min)
hansYr = hansYr+1 # Hansen Year starts at 0 so I want to add 1 to all
colnames(hansYr) = paste0(colnames(hansYr), "_30")
# file = paste0("data2/Prediction_150mGrid_GRNP_HansenYear_maxMin_30m.csv")
# file = paste0("data2/Prediction_150mGrid_buffer_HansenYear_maxMin_30m.csv")
# write.csv(x = hansYr, file = file, row.names = F)


# 3.Calculate Buffer Values for the 3 GLCM statistics   
p = spTransform(spdf, crs(sat))
lsat_mean = extract(x = sat,
                    y = p,
                    buffer=30,
                    fun=mean,
                    df=F)
colnames(lsat_mean) = paste0(colnames(lsat_mean), "_30")

lsat_SD = extract(x = sat,
                  y = p,
                  buffer=30,
                  fun=sd,
                  df=F)
colnames(lsat_SD) = paste0(colnames(lsat_SD), "_30")
glcm_df = cbind(lsat_mean, lsat_SD)
# filename2 = paste0("data2/Prediction_150mGrid_GRNP_GLCM_Landsat_2022_Jan_30m.csv")
# filename2 = paste0("data2/Prediction_150mGrid_buffer_GLCM_Landsat_2022_Jan_30m.csv")
# write.csv(glcm_df, filename2, row.names = F)  


#############################################
# These data are scaled in the model & therefore -
# As per Li - Clamp the Values to the Max & Min values at the trap sites Before Making Canonical Components

# Radar values
tru_rad = read.csv("data2/golaTRAPpoints_Sentinel1_2022_q1_30_500m.csv")
tru_rad = tru_rad[grep("30", x = colnames(tru_rad))]
max = apply(tru_rad, MARGIN = 2, FUN = max)
min = apply(tru_rad, MARGIN = 2, FUN = min)

for (i in 1:ncol(radar)) {
    radar[,i] = ifelse(radar[,i] > max[[i]],max[[i]],radar[,i])
    radar[,i] = ifelse(radar[,i] < min[[i]],min[[i]],radar[,i])
}

# GLCM Values
tru_glcm = read.csv("data2/golaTRAPpoints_GLCM_Landsat_2022_Jan_30-500m.csv")
tru_glcm = tru_glcm[grep("30", x = colnames(tru_glcm))]
max = apply(tru_glcm, MARGIN = 2, FUN = max)
min = apply(tru_glcm, MARGIN = 2, FUN = min)

for (i in 1:ncol(glcm_df)) {
  glcm_df[,i] = ifelse(glcm_df[,i] > max[[i]],max[[i]],glcm_df[,i])
  glcm_df[,i] = ifelse(glcm_df[,i] < min[[i]],min[[i]],glcm_df[,i])
}

#####################################################################
## Use these variables to Create Canonical Components
library(sgdm)
library(gdm)
library(dplyr)

## Load the vectors to create the Canonical Components from
rad_vect = read.csv("data2/sgdm/vectors_to_create_components/vectors_for_Can_Cpnts_k1_RADAR_30.csv", row.names = 1)
lsat_vect = read.csv("data2/sgdm/vectors_to_create_components/vectors_for_Can_Cpnts_k10_lsat_glcm-bandslsat_glcm-indices_30.csv", row.names = 1)

## Radar ##
# Put band names in order for both vector & values 
radar = radar[,order(colnames(radar))]
rad_vect = rad_vect[order(rownames(rad_vect)),]
rad_vect = as.matrix(rad_vect[,-2])

# Create a df to add the components into
rad_can = matrix(data = NA, nrow = nrow(radar), ncol = 1) # Columns = Number of Components
rad_can = as.data.frame(rad_can)

# Multiply each counted value by its vector value, then sum them to create the single canonical component for radar
for (i in 1:nrow(rad_can)) {
  for (j in 1:1) {
    rad_can[i,j] = sum(rad_vect[,j]*radar[i,])
  }
}
names(rad_can)[1] = "rad_can_30"
# write.csv(rad_can, "data2/sgdm/components/Prediction_150mGrid_GRNP_Rad_K1_30m.csv")
# write.csv(rad_can, "data2/sgdm/components/Prediction_150mGrid_buffer_Rad_K1_30m.csv")

## Satellite ##
# Put band names in order for both vector & values 
glcm_df = glcm_df[,order(colnames(glcm_df))]
lsat_vect = lsat_vect[order(rownames(lsat_vect)),]
lsat_vect = lsat_vect[,-11]

# Create a df to add the components into
lsat_can = matrix(data = NA, nrow = nrow(glcm_df), ncol = 10) # Columns = Number of Components
lsat_can = as.data.frame(lsat_can)

# Multiply each counted value by its vector value, then sum them to create the single canonical component for radar
for (i in 1:nrow(lsat_can)) {
  for (j in 1:ncol(lsat_can)) {
    lsat_can[i,j] = sum(lsat_vect[,j]*glcm_df[i,])
  }
}
names(lsat_can) = c("cc1_3-",'cc2_30','cc3_30','cc4_30','cc5_30','cc6_30','cc7_30','cc8_30','cc9_30','cc1_30')
# write.csv(lsat_can, "data2/sgdm/components/Prediction_150mGrid_GRNP_glcmB_glcmI_K10_30m.csv")
# write.csv(lsat_can, "data2/sgdm/components/Prediction_150mGrid_buffer_glcmB_glcmI_K10_30m.csv")



