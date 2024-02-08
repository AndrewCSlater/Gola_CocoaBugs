library(recluster)
library(Rtsne)
library(ggplot2)
library(sf)
######### Load Names of Best Predicted Species from Cross Validation of Tuned Model
# pred.mn = read.csv("Data/Predicted_Occurrence_150mGrid.csv")
pred.mn = read.csv("data2/Predicted_Occurrence_150mGrid_BufferGRNP.csv")
best_sp = read.csv("data2/Results_AUC_perSp_Final_Model.csv")
use_sp = best_sp[which(best_sp$times.hi.auc>=2) ,1]

XY_in = read.csv("data2/prediction_coordinate_150mgrid_inside_GRNP.csv")
XY_in = XY_in[,(4:6)]
XY_buf = read.csv("data2/prediction_coordinate_150mgrid_gola_buffer.csv")
XY_buf = XY_buf[,c(1,4:5)]
names(XY_in) = c("SITE_ID", "Long", "Lat")
names(XY_buf) = c("SITE_ID", "Long", "Lat")
XY = rbind(XY_in, XY_buf)
xy = XY[,2:3]


incidence = best_sp[which(best_sp$times.hi.auc>=2) ,13]

to_map = as.data.frame(pred.mn[,colnames(pred.mn)%in%use_sp])

Xmat = to_map

perplexity = 140			# Rule of thumb = (Number of Sites)^0.5


# ## Initial PCA dimensions
# # https://towardsdatascience.com/how-to-tune-hyperparameters-of-tsne-7c0596a18868

tsne = Rtsne::Rtsne(Xmat, dims = 2, perplexity = perplexity, theta = 0.5, pca = FALSE, num_threads = 0) # can set parallel options if using openMP

col <- recluster::recluster.col(tsne$Y) # Create & add RGB values based on 2 PCA values

leakage = st_read("data2/shapefiles/Leakage_Belt_clipped_to_1kmbufferEnvelope_for_Predictions.shp")
leakage = st_transform(leakage, CRS("+proj=longlat +datum=WGS84 +no_defs"))

survey = read.csv("data2/XY_No_Satellite_Gola.csv", row.names = 1)
rownames(survey) = 1:nrow(survey)

recluster.plot.col(mat = col)
recluster.plot.sites.col(long = xy$Lon, lat = xy$Lat, mat = col, text = F, pch = 20, cex=0.8,  cex.main = 0.85, main = "",  ylab = "", xlab = "") # Plot the tsne RGB values
points(survey[,c(6,5)], pch=19, cex=0.8) # Add the trap points for reference
plot(leakage$geometry, lwd=3, add=T)
# Add North Arrow & Scale bar
prettymapr::addnortharrow(pos = "bottomright", scale = 0.5, cols = "darkgrey", padin = c(0.4, 0.2))
prettymapr::addscalebar(plotunit="latlon", plotepsg = 2162, pos = "topleft", htin = 0.05, bar.cols = c("darkgrey", "white"))
