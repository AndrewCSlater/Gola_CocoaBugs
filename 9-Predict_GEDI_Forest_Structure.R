# rgdal: for vector work; sp package should always load with rgdal. 
library(rgdal)
# raster: for metadata/attributes- vectors or rasters
library (raster)
library(sf)
library(dplyr)
library(ggplot2)

#### CREATE MODEL TO PREDICT FROM ####

## Load Data
gedi = read.csv("data2/Gedi_in_Gola.csv", row.names = 1)
lsat = read.csv("data2/sgdm/components/GEDI_Prediction_golaBufferZone_glcmB_glcmI_K10_30m.csv", row.names = 1)
rad = read.csv("data2/sgdm/components/GEDI_Prediction_golaBufferZone_Rad_K1_30m.csv", row.names = 1)
hans = read.csv("data2/GEDI_prediction_golaBufferZone_HansenYear_30m.csv")
hans = hans[,1]
dist = read.csv("data2/GEDI_prediction_golaBufferZone_Distance_to_GRNP_boundary.csv")
dist = dist[,4]

# Create DF that will be used to Train Model & Test Predictions
df_pred = cbind(lsat,hans,dist,rad)

####################################
### --- DO RANDOM FOREST etc --- ###
####################################

# Turn points into a Spatial Object with a Generic CRS

xy = gedi[,1:2]
spdf <- SpatialPointsDataFrame(coords = xy, data = gedi)


# Initiate Random Forest for analysis
library(randomForest)
library(caret)
library(caTools)

# Create df with:
# - Response Var in Column Height=4, Pct Canopy Cover=9, Height Diversity=6
# - Only Predictors being used in furthers columns

# There are some -9999 Diversity & Cover values which means the readings are suspect, & so I remove those rows from analysis
gedi = gedi[,c(1,2,4,6,9)]
df_mod = cbind(gedi,df_pred)
df_mod = df_mod[df_mod$Height_Diversity>=0 & df_mod$canopy_cover_pct>=0,] ### Excludes -9999 values

# Split data into Training & Testing
sample <- sample.split(df_mod[,1], SplitRatio = 0.75)
train  <- subset(df_mod, sample == TRUE)
test   <- subset(df_mod, sample == FALSE)

## Create Train & Test df for each of the 3 structures
train_height = train[,c(3,6:ncol(train))]
train_height = train_height[!is.na(train_height[,1]),]
test_height = test[,c(3,6:ncol(test))]

train_height_div = df_mod[,c(4,6:ncol(df_mod))]
train_height_div = train_height_div[!is.na(train_height_div[,1]),]
test_height_div = test[,c(4,6:ncol(test))]

train_cover = df_mod[,c(5,6:ncol(df_mod))]
train_cover = train_cover[!is.na(train_cover[,1]),]
test_cover = test[,c(5,6:ncol(test))]

### Run RF models for each structure
mod_height = randomForest(train_height$height_m ~., data = train_height, mtry = 4, ntree = 2001, importance=T)
mod_height_div = randomForest(train_height_div$Height_Diversity ~., data = train_height_div, mtry = 4, ntree = 2001, importance=T)
mod_cover = randomForest(train_cover$canopy_cover_pct ~., data = train_cover, mtry = 4, ntree = 2001, importance=T)

### Plot variables by importance to model
varImpPlot(mod_height)
varImpPlot(mod_height_div)
varImpPlot(mod_cover)

### Predict each model structure using test data
result_height = data.frame(test_height$height_m, predict(mod_height, test_height[,2:ncol(test_height)], type="response"))
names(result_height) = c("Measured_Height", "Predicted_Height")
result_height_div = data.frame(test_height_div$Height_Diversity, predict(mod_height_div, test_height_div[,2:ncol(test_height_div)], type="response"))
names(result_height_div) = c("Measured_Height-Div", "Predicted_Height-Div")
result_cover = data.frame(test_cover$canopy_cover_pct, predict(mod_cover, test_cover[,2:ncol(test_cover)], type="response"))
names(result_cover) = c("Measured_Cover", "Predicted_Cover")

### Plot Predicted by Measured structure values
# Linear models of Predicted by Measured structure
lin_height = lm(result_height$Predicted_Height~result_height$Measured_Height)
summary(lin_height)
plot(result_height)
abline(lin_height)

lin_height_div = lm(result_height_div$`Predicted_Height-Div`~result_height_div$`Measured_Height-Div`)
summary(lin_height_div)
plot(result_height_div)
abline(lin_height_div)

lin_cover = lm(result_cover$Predicted_Cover~result_cover$Measured_Cover)
summary(lin_cover)
plot(result_cover)
abline(lin_cover)
