library(dplyr)
library(raster)
library(sf)
library(rgdal)

####################################
# --- Load the Trap Point Data
site = read.csv("data2/X_Site_Env_data.csv")
traps = site[,c(3, 6,5)]

xy = traps[,2:3]
spdf <- SpatialPointsDataFrame(coords = xy, data = traps,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

# # ----------------------------
# # Import the satellite images
# # ----------------------------
rad = brick("satellite_data/Sent1_2022_q1_Radar.tif")
hansen = raster("satellite_data/Hansen_GFC-2019-v1.7_lossyear_10N_020W.tif")
lsat = stack("satellite_data/Landsat8_2022_Jan_Bands_Indices_terra.tif")

b = c(30, 60, 90, 120, 150, 250, 500) # Buffer radii

### CREATE VALUES FOR EACH BUFFER RADIUS ###
# 1. Radar - Sentinel-1 
radar = c()
for (i in b) {
# i=30
# Sentinel_1 VV & VH Radar
p = spTransform(spdf, crs(rad))
rad_mean = extract(x = rad,
                     y = p,
                     buffer=i,
                     fun=mean,
                     df=F)
  
rad_SD = extract(x = rad,
                   y = p,
                   buffer=i,
                   fun=sd,
                   df=F)
s1 = cbind(rad_mean, rad_SD)
colnames(s1) = c("VV_mean", "VH_mean", "VV_sd", "VH_sd")
colnames(s1) = paste0(colnames(s1), "_", i)
radar = cbind(radar, s1)
}
# file = paste0("data2/golaTRAPpoints_Sentinel1_2022_q1_30_500m.csv")
# write.csv(x = radar, file = file, row.names = F)

# 2. Hansen year of deforesetation - NB I use Max & Min of year within the buffer
hansYr = c()
for (i in b) {
# Hansen Year of Last Logging
p = spTransform(spdf, crs(hansen))
hansYr_max = extract(x = hansen,
                       y = p,
                       buffer = i,
                       fun = max,
                       df=F)
  
hansYr_min = extract(x = hansen,
                       y = p,
                       buffer = i,
                       fun = min,
                       df=F)
h = cbind(hansYr_max, hansYr_min)
h = h+1 # Hansen Year starts at 0 so I want to add 1 to all
colnames(h) = paste0(colnames(h), "_", i)
hansYr = cbind(hansYr, h)
}
# file = paste0("data2/golaTRAPpoints_HansenYear_maxMin_30_500m.csv")
# write.csv(x = hansYr, file = file, row.names = F)

# 3. Landsat
lst8 = c()
for (i in b) {
# Landsat Satellite Imagery
p = spTransform(spdf, crs(lsat))
lsat_mean = extract(x = lsat,
                      y = p,
                      buffer=i,
                      fun=mean,
                      df=F)
colnames(lsat_mean) = paste0(colnames(lsat_mean),"_mean_",i)
lsat_SD = extract(x = lsat,
                    y = p,
                    buffer=i,
                    fun=sd,
                    df=F)
colnames(lsat_SD) = paste0(colnames(lsat_SD),"_SD_",i)
  
lst8 = cbind(lst8, lsat_mean, lsat_SD)
}
# file = paste0("data2/golaTRAPpoints_Landsat_2022_Jan_30_500m.csv")
# write.csv(x = lst8, file = file, row.names = F)

##########################################################
### CREATE GLCM VALUES FOR ALL BANDS & INDICES & BUFFER ###
library(glcm)

### Creating GLCM over a huge raster may be excessively resource heavy, so crop the rasters to the bounding envelope of the buffered points.
# import the Gola Buffer Envelope vector boundary
crop = readOGR("data2/shapefiles/Gola_MalaiseTrap_Buff1000Envelope_.shp")

crop = spTransform(crop, crs(lsat))
lsat_crop = crop(lsat, crop)

glcm_df=c()

# 1. Create a single layer raster from of each band of the Lsat image in turn 
for (band in c(1:length(names(lsat_crop)))) { # Chosen bands in the Raster Brick
  ST = Sys.time()
  img = lsat_crop[[band]]
  n = names(img)
  
# 2.Create a Raster stack of the 3 GLCM statistic values based on the current raster band
glcm = glcm(img,
            window = c(3,3), 
            shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
            n_grey = 16, # 16 = 4-Bit, 64 = 6-Bit, 256 = 8-Bit
            statistics = c("mean", "contrast", "entropy"))

# 3.Calculate Buffer Values for the 3 GLCM statistics   
p = spTransform(spdf, crs(glcm))
for (j in b) {
lsat_mean = extract(x = glcm,
                    y = p,
                    buffer=j,
                    fun=mean,
                    df=F)
colnames(lsat_mean) = paste0("lsat_", n,"_mean_",j,"m_", colnames(lsat_mean))
  
lsat_SD = extract(x = glcm,
                  y = p,
                  buffer=j,
                  fun=sd,
                  df=F)
colnames(lsat_SD) = paste0("lsat_", n,"_SD_",j,"m_", colnames(lsat_SD))
  
glcm_df = cbind(glcm_df, lsat_mean, lsat_SD)
} # End Buffer extraction of GLCM

# filename2 = paste0("data2/golaTRAPpoints_GLCM_Landsat_2022_Jan_30-500m.csv")
# write.csv(glcm_df, filename2, row.names = F)  

} # End creation of all GLCM

####################################################################
### CREATE CANONICAL COMPONENTS FOR ALL GROUPS and BUFFER LEVELS ###
library(sgdm)
library(gdm)
library(dplyr)

# Load Survey Data
survey = read.csv("data2/XY_No_Satellite_Gola.csv", row.names = 1)
rownames(survey) = 1:nrow(survey)
survey$SITE_ID = rownames(survey)

# Create Count Data
y = survey[,c(20:ncol(survey))]
y[y>0]=1

# As per Pichler 2021 - Remove OTU's with fewer than 3 occurrences
# He Removes Sites with fewer than 4 OTU's - But I keep all sites that have any insects
df1 = survey
r = which(rowSums(y)>0) # Rows to keep
df1 = df1[c(r),]
c = which(colSums(y)>2) # Columns to keep
c = c+19
df1 = df1[,c(1:19,c)]
y = df1[,c(20:ncol(df1))]
y[y>0]=1

# Load Satellite data
lsat = read.csv("data2/golaTRAPpoints_Landsat_2022_Jan_30_500m.csv")
lsat = lsat[rownames(lsat)%in%rownames(df1),]

sat_glcm = read.csv("data2/golaTRAPpoints_GLCM_Landsat_2022_Jan_30-500m.csv")
sat_glcm = sat_glcm[rownames(sat_glcm)%in%rownames(df1),]

radar = read.csv("data2/golaTRAPpoints_Sentinel1_2022_q1_30_500m.csv")
radar = radar[rownames(radar)%in%rownames(df1),]

# Create Predictor (bands) and Response (count) Data Sets
df2 = df1[, c("SITE_ID", "Longitude", "Latitude")]
df2 = rename(df2, Plot_ID = SITE_ID, X=Longitude, Y=Latitude)

df2$Plot_ID = as.numeric(df2$Plot_ID) # NB: - The plot_id needs to be numeric and unique

insect = cbind(df2$Plot_ID, y) # Include plot_id in 1st column
colnames(insect)[1] = "Plot_ID"

# Split data into Bands & Indices
band_names = c("blue", "green", "red", "NIR", "SWIR", "MIR", "RE1", "RE2", "RE3")

b = c(30, 60, 90, 120, 150, 250, 500) # Buffer radii

for (buffer in b) {
#buffer=30
# Filter by buffer size
raw_buffer = lsat[,grep(paste(buffer,collapse="|"),colnames(lsat))]
glcm_buffer  = sat_glcm[,grep(paste(buffer,collapse="|"),colnames(sat_glcm))]
###   OR
# raw_buffer = radar[,grep(paste(buffer,collapse="|"),colnames(radar))]

# Divide Raw & GLCM by Band & Index
raw_bands = raw_buffer[,grep(paste(band_names,collapse="|"),colnames(raw_buffer))]
raw_indices = raw_buffer[,-grep(paste(band_names,collapse="|"),colnames(raw_buffer))]
glcm_bands = glcm_buffer[,grep(paste(band_names,collapse="|"),colnames(glcm_buffer))]
glcm_indices = glcm_buffer[,-grep(paste(band_names,collapse="|"),colnames(glcm_buffer))]

satellite_used = list(raw_bands, raw_indices, glcm_bands, glcm_indices)
satellite_used_names = c("lsat_bands", "lsat_indices", "lsat_glcm-bands", "lsat_glcm-indices")

# Create all combinations of groups - I already have individual Bands & Indices & combined glcm 0f the 15 possible combinations, so am creating the other 12
a = list(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]], satellite_used[[4]], cbind(satellite_used[[1]], satellite_used[[2]]), cbind(satellite_used[[1]], satellite_used[[3]]), cbind(satellite_used[[1]], satellite_used[[4]]), cbind(satellite_used[[2]], satellite_used[[3]]), cbind(satellite_used[[2]], satellite_used[[4]]), cbind(satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[2]], satellite_used[[3]], satellite_used[[4]]), cbind(satellite_used[[1]], satellite_used[[2]], satellite_used[[3]], satellite_used[[4]]))

bb = c(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]], paste0(satellite_used_names[[1]], satellite_used_names[[2]]), paste0(satellite_used_names[[1]], satellite_used_names[[3]]), paste0(satellite_used_names[[1]], satellite_used_names[[4]]), paste0(satellite_used_names[[2]], satellite_used_names[[3]]), paste0(satellite_used_names[[2]], satellite_used_names[[4]]), paste0(satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]]), paste0(satellite_used_names[[1]], satellite_used_names[[2]], satellite_used_names[[3]], satellite_used_names[[4]]))


for (s in 1:15) {  ## For each combination of satellite data
#  s=15
sat = a[[s]]

# There must not be any NA in either matrix
# This replaces any NA in the satellite data with the column (Band) mean value
for(ii in 1:ncol(sat)){
  sat[is.na(sat[,ii]), ii] <- mean(sat[,ii], na.rm = TRUE)
}

bbb = bb[[s]] # Current Group Name


########################################
# Create 1 Canon from the 4 Radar values
# sat = raw_buffer
  
# The band data must have variation and cannot have a zero Standard Deviation
# Find columns with SD of 0 (There is no variability in the column)
r = which(apply(sat,2,sd)==0)

use = if(is.na(r[2])) {
  sat
} else {
  sat[,-c(r)]
}

bands = cbind(df2, use) # Include plot_id and X Y in first 3 columns

### RUN THE SGDM PROCESS
# 1 - Test models using a range of penalisation values from 0.6 - 1.0 on both 
# predictor and response variables - the output is a 5x5 vector of RMSE values
mod = sgdm.param(predData = bands, bioData = insect, k = 10)
###   OR
# mod = sgdm.param(predData = bands, bioData = insect, k = 1) ### FOR RADAR
gc()
# file = paste0("data2/sgdm/ModParams_gola_k10_",bbb,"_",buffer,".csv")
# file = paste0("data2/sgdm/ModParams_gola_k1_RADAR_",buffer,".csv")
# write.csv(mod, file = file)

# 2 - Retrieve the best performing model based on the highest vector value -
# (using the 'sgdm.best' function, with  output option 'm')
#  This provides the appropriate penalisation values to be used in the model
best_mod = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "m", k = 10)
# best_mod = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "m", k = 1) ### FOR RADAR

# 3 - Extract Sparse Canonical Components relating to the best model
# (using the 'sgdm.best' function, with  output option 'c')
best_mod_cpnts = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "c", k = 10)
best_mod_cpnts = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "c", k = 1) ### FOR RADAR
# file = paste0("data2/sgdm/components/Can_Cpnts_k10_",bbb,"_",buffer,".csv")
# file = paste0("data2/sgdm/components/Can_Cpnts_k1_RADAR_",buffer,".csv")
# write.csv(best_mod_cpnts, file = file, row.names = F)

## Check variable contributions
# 8 - Extract the canonical vectors of the best model -
# (using the 'sgdm.best' function, with output option 'v')
best_mod_vects = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "v", k = 10)
best_mod_vects = sgdm.best(perf.matrix = mod, predData = bands, bioData = insect, output = "v", k = 1) ### FOR RADAR

# Add the band names to the vector
best_mod_vects = as.data.frame(best_mod_vects)
rownames(best_mod_vects) = colnames(bands[,-c(1:3)])
best_mod_vects$band = colnames(bands[,-c(1:3)])

# Save this value to use to calculate Canonical Components across all pixels of an image
# file = paste0("data2/sgdm/vectors_to_create_Can_Cpnts_k10_",bbb,"_",buffer,".csv")
# file = paste0("data2/sgdm/vectors_to_create_Can_Cpnts_k1_RADAR_",buffer,".csv")
# write.csv(best_mod_vects, file = file, row.names = T)
gc()
} # End of Each Satellite Group Loop
} # End of Each Buffer Size 



