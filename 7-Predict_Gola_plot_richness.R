# rgdal: for vector work; sp package should always load with rgdal. 
library(rgdal)
# raster: for metadata/attributes- vectors or rasters
library (raster)
library(sf)
library(sjSDM)
library(dplyr)
library(ggplot2)


#### 1 - LOAD MODEL TO PREDICT FROM ####
load("data2/models/Final_K10.rad.hans.dist.Rdata")


### Load Data
## Survey Data - Counts & Environment
survey = read.csv("data2/XY_No_Satellite_Gola.csv", row.names = 1)
rownames(survey) = 1:nrow(survey)

## Data to Predict From 
XY_in = read.csv("data2/prediction_coordinate_150mgrid_inside_GRNP.csv")
XY_in = XY_in[,(4:6)]
XY_buf = read.csv("data2/prediction_coordinate_150mgrid_gola_buffer.csv")
XY_buf = XY_buf[,c(1,4:5)]
names(XY_in) = c("SITE_ID", "Long", "Lat")
names(XY_buf) = c("SITE_ID", "Long", "Lat")
XY = rbind(XY_in, XY_buf)
coords = XY

lsat_in = read.csv("data2/sgdm/components/Prediction_150mGrid_GRNP_glcmB_glcmI_K10_30m.csv", row.names = 1)
lsat_buf = read.csv("data2/sgdm/components/Prediction_150mGrid_buffer_glcmB_glcmI_K10_30m.csv", row.names = 1)
lsat = rbind(lsat_in, lsat_buf)
rad_in = read.csv("data2/sgdm/components/Prediction_150mGrid_GRNP_Rad_K1_30m.csv", row.names = 1)
rad_buf = read.csv("data2/sgdm/components/Prediction_150mGrid_buffer_Rad_K1_30m.csv", row.names = 1)
rad = rbind(rad_in, rad_buf)
hans_in = read.csv("data2/Prediction_150mGrid_GRNP_HansenYear_maxMin_30m.csv")
hans_buf = read.csv("data2/Prediction_150mGrid_buffer_HansenYear_maxMin_30m.csv")
hans = rbind(hans_in, hans_buf)
hans = hans[,1]
dist_in = read.csv("data2/prediction_150mGrid_GRNP_DistancetoBoundary_from_INSIDE.csv")
dist_buf = read.csv("data2/Prediction_150mGrid_buffer_DistanceToPark.csv")
dist = c(dist_in[,4]*-1, dist_buf[,4])

# Create DF that will be used to Fit Model
df_pred = cbind(lsat,hans,dist,rad)
names(df_pred) = model.train$names[-1]
# There must not be any NA, this replaces NA with the column mean value
for(i in 1:ncol(df_pred)){
  df_pred[is.na(df_pred[,i]), i] <- mean(df_pred[,i], na.rm = TRUE)
}

# Spatial terms: for linear spatial model
# XY_pred = coords[,4:5]
XY_pred = coords[,2:3]
names(XY_pred) = c("X1", "X2")

#### Scale Data
newdd = scale(df_pred)
newsp = scale(XY_pred)

#############################################################
### MAKE PREDICTIONS ###
## predict 5 times and save as list
pred_all <- lapply(1:5, function(i) {
  predict(model.train, newdata = newdd, SP = newsp)}
)

# get mean prediction
pred.mn = simplify2array(pred_all)
pred.mn = apply(pred.mn, MARGIN = c(1,2), FUN = mean, na.rm=T)
colnames(pred.mn) = model.train$species
# write.csv(pred.mn, "data2/Predicted_Occurrence_150mGrid_BufferGRNP.csv", row.names = F)

##############################################################

######### Load Names of Best Predicted Species from Cross Validation of Final Model
pred.mn = read.csv("data2/Predicted_Occurrence_150mGrid_BufferGRNP.csv")
best_sp = read.csv("data2/Results_AUC_perSp_Final_Model.csv")
use_sp = best_sp[which(best_sp$times.hi.auc>=2),1] # Select those species with a hi-auc in at least 2 of the runs

xy = coords[,2:3]

incidence = best_sp[which(best_sp$times.hi.auc>=2) ,13]
to_map = as.data.frame(pred.mn[,colnames(pred.mn)%in%use_sp]) # Select only species with high AUC to plot

#################################################
## To Plot Richness (Of Best Predicted Species)

### Sum probabilities of occurrence of all species at each predicted site to get its predicted richness
to_map$rich = rowSums(to_map)
to_map = cbind(xy,to_map)

leakage = st_read("data2/Shapefiles/Leakage_Belt_clipped_to_1kmbufferEnvelope_for_Predictions.shp")
leakage = st_transform(leakage, CRS("+proj=longlat +datum=WGS84 +no_defs"))

gp = ggplot(to_map, aes(Long, Lat))

gp2 = gp + geom_point(aes(colour=rich), size = 0.5) +
  scale_colour_gradient2(low = "darkblue", mid = "gold", high = "green", midpoint = 15) +
  labs(colour = "Species\nRichness") +
  theme(panel.background = element_rect(fill='snow'))

gp2 + geom_sf(data = leakage, inherit.aes = FALSE, fill = NA, color = "black", lwd=1.25) +
  geom_point(data = survey, mapping = aes(Longitude, Latitude),  inherit.aes = FALSE, pch=19, cex=0.8) +
# Add scale and North arrow
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("darkgrey", "white"),
    text_family = "ArcherPro Book"
  ) +
  ggspatial::annotation_north_arrow(
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    location = "br", which_north = "true",
    pad_x = unit(0.5, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_orienteering(
      fill = c("darkgrey", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"))

