library(readxl)
library(dplyr)

####################################
# --- Read in Site Data
site = readxl::read_excel(path = "data2/CocoaBugs_field_data_20220813.xlsx", sheet = "Section_all") 

str(site)
# Entry Order not required
site = within(data = site, expr = rm(Entry_order))

# Surveyor, Site ID, Stumps, Burn, Oil Palm need to be as Factor
fac = c("Surveyor", "recent_burning_y_n", "stumps_gt_10cm__y_n", "stumps_gt_30cm_y_n", "oil_palm_y_n")
site[fac] <- lapply(site[fac], as.factor)
rm(fac)

# Edit Site ID to be the same as for the count data
site$SITE_ID = ifelse(site$Section=="GRNP", paste0("GRNP",site$SITE_ID), site$SITE_ID)
site$SITE_ID = gsub(pattern = " ", replacement = "_", x = site$SITE_ID)
site$SITE_ID = gsub(pattern = "_new", replacement = "", x = site$SITE_ID)

site = within(data = site, expr = rm(Section))

# Time deployed & collected are not required
# All deployments are for 5 days apart from 1 of 4 days range(site$date_collect-site$date_deploy) -
# - and all over a 6 week period range(site$date_deploy) - so dates are unlikely to be of any use
site = within(data = site, expr = rm(date_deploy, time_deploy, datetime_deploy, date_collect, time_collect, datetime_collect))

# Lat & Lon need to be Numeric
site$Latitude = as.numeric(site$Latitude)
site$Longitude = as.numeric(site$Longitude)
# NA have been identified on 1 site - site had to be moved and coordinates supplied independently
# Install them here 
s67 = readxl::read_excel(path = "data2/CocoaBugs_field_data_20220813.xlsx", sheet = "Gola_selected_sites_WGS1984")
s67Lat = s67$Latitude[s67$SITE_ID==67]
s67Lon = s67$Longitude[s67$SITE_ID==67]
site$Latitude[site$SITE_ID=="67 new"] = s67Lat
site$Longitude[site$SITE_ID=="67 new"] = s67Lon
rm(s67, s67Lat, s67Lon)
gc()

# Canopy ht, bare ground, grass, cocoa, liana, forest, agri need to be converted to numeric
names(site)
site[,5:11] = apply(site[,5:11], 2, function(x) gsub(pattern = " ", replacement = "", x))
site[,5:11] = apply(site[,5:11], 2, function(x) gsub(pattern = "%", replacement = "", x))
site[,5:11] = apply(site[,5:11], 2, function(x) gsub(pattern = ">", replacement = "", x))
site[,5:11] = apply(site[,5:11], 2, function(x) gsub(pattern = "<", replacement = "", x))
site[,5:11] = apply(site[,5:11], 2, function(x) gsub(pattern = "m", replacement = "", x))

# site[,5:11] <- lapply(site[,5:11], as.numeric)
site[,5:11] <- lapply(site[,5:11], as.factor)

# Audiomoth & BAR recorder not relevant unless we get sound data
site = within(data = site, expr = rm(Audiomoth, BAR_Recorder))

# Delete rows with NO sample tube & then delete sample tube column
site = site[site$sample_tube=="present",]
site = within(data = site, expr = rm(sample_tube))

# write.csv(site, "data2/X_Site_Env_data.csv")

###############################
# --- Read in Survey Data

### ID - Link between Count & Site codes
id = readxl::read_excel(path = "data2/Malaise_traps_ZG0001592.client.otuids.xlsx", sheet = "SO00140.Leray2.QC")

### Individual Taxonomic Unit ("species") counts
sp = readxl::read_excel(path = "data2/Malaise_traps_ZG0001592.client.otuids.xlsx", sheet = "ZG0001592.client.otuids")

str(sp)
sp = sp[!is.na(sp$Phylum),]        # Remove rows with NA ID
sp = sp[sp$Class=="Insecta",]       # Select only ITU identified as Insecta
sp = sp[!is.na(sp$Class),]
sp1 = as.data.frame(t(sp[, 19:length(sp)])) # Transpose only count data
rownames(sp1) = gsub("_MiSe.*", "", rownames(sp1)) # Reduce Row Names to equate to DNA ID from the
rownames(sp1) = gsub(".*_", "", rownames(sp1))     # ID sheet to then get Site Name for counts

sp1[is.na(sp1)]=0 # There are some NA values in the counts, which I change to value 0
sp1$DNA_ID = rownames(sp1)
sp1 = select(sp1, ncol(sp1), everything())

str(id)
# Columns that are not required
id = within(data = id, expr = rm(Library_type, Date_received, Kit_ID, Volume_filtered,Qubit,MiSeq,PCR_notes,Latitude, Longitude, Marker, Company, Project_name, Sample_type))
id = id[,-5]
colnames(id) = gsub(" ", "_", colnames(id))


############################
# Add Counts to Sites
y = right_join(id, sp1, by="DNA_ID")

y1 =y
y1$DNA_ID = gsub("[^0-9.-]", "", y1$DNA_ID) # Remove the split DNA ID to group & Sum them
y2 = aggregate(y1[,5:ncol(y1)], y1[,1:2], FUN = sum ) # Groups unique sets of DNA_ID and Site name and sums the sightings for all species

# Edit the Site Codes column
y2$Customer_ID = gsub(pattern = "Section ", replacement = "", x = y2$Customer_ID)
y2$Customer_ID = gsub(pattern = " NEW", replacement = "_new", x = y2$Customer_ID)
y2$Customer_ID = gsub(pattern = "GRNP-", replacement = "GRNP", x = y2$Customer_ID)
y2$Customer_ID = gsub(pattern = " ", replacement = "_", x = y2$Customer_ID)
y2$Customer_ID = gsub(pattern = "-", replacement = "_", x = y2$Customer_ID)
names(y2)[2] = "SITE_ID"
y2$SITE_ID = gsub(pattern = "_new", replacement = "", x = y2$SITE_ID)
y2$SITE_ID = gsub("^[^_]*_","", y2$SITE_ID) # Remove everything before & including the first underscore - matching zero or more characters that are not an underscore ([^_]*) from the start (^) of the string, followed by an underscore (_)

# Sites named 66,67 and 68,81 are just named 66 and 68 in X Data so change the names here to comply with that
y2$SITE_ID = gsub(pattern = ",67", replacement = "", x = y2$SITE_ID)
y2$SITE_ID = gsub(pattern = ",81", replacement = "", x = y2$SITE_ID)
y2$SITE_ID = gsub(pattern = "Waypoint_", replacement = "", x = y2$SITE_ID)

# write.csv(y2, "data2/Y_Counts_per_site.csv", row.names = F)

####
# - Bind the Environment & Count data into 1 DF for import into analysis or to have satellite data added

df = right_join(site, y2, by = "SITE_ID")

# Put Environmental data into numerical or ordered factoral structure and turn percentages into decimal
names(df)
env = df[,c(2:17)]
# str(env)
# fac = c("SITE_ID")
# env[fac] <- lapply(env[fac], as.factor)

env$bare_ground_pct = as.numeric(env$bare_ground_pct)
env$bare_ground_pct = ifelse(env$bare_ground_pct >1, env$bare_ground_pct/100, env$bare_ground_pct)
env$bare_ground_pct[is.na(env$bare_ground_pct)] = 0

env$agric_pct
env$agric_pct = ifelse(env$agric_pct==">60%",60,env$agric_pct)
env$agric_pct = as.numeric(env$agric_pct)
env$agric_pct = ifelse(env$agric_pct >1, env$agric_pct/100, env$agric_pct)

env$oil_palm_y_n = ifelse(env$oil_palm_y_n == "Yes",1,0)

env$stumps_gt_30cm_y_n = ifelse(env$stumps_gt_30cm_y_n == "Yes",1,0)
env$stumps_gt_10cm__y_n = ifelse(env$stumps_gt_10cm__y_n == "Yes",1,0)

env$recent_burning_y_n = ifelse(env$recent_burning_y_n == "Yes",1,0)

unique(env$forest_pct)
env$forest_pct = ifelse(env$forest_pct == "5-20",12.5,env$forest_pct)
env$forest_pct = ifelse(env$forest_pct == "May-20",12.5,env$forest_pct)
env$forest_pct = ifelse(env$forest_pct=="20-40",30,env$forest_pct)
env$forest_pct = ifelse(env$forest_pct=="40-60",50,env$forest_pct)
env$forest_pct = as.numeric(env$forest_pct)
env$forest_pct = ifelse(env$forest_pct >1, env$forest_pct/100, env$forest_pct)

unique(env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "5-10",7.5,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "05-Oct",7.5,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "5-20",12.5,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "May-20",12.5,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "10-20",15,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct == "Oct-20",15,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct=="20-40",30,env$liana_pct)
env$liana_pct = ifelse(env$liana_pct=="40-60",50,env$liana_pct)
env$liana_pct = as.numeric(env$liana_pct)
env$liana_pct = ifelse(env$liana_pct >1, env$liana_pct/100, env$liana_pct)

unique(env$cocoa_pct)
env$cocoa_pct = ifelse(env$cocoa_pct == "5-20",12.5,env$cocoa_pct)
env$cocoa_pct = ifelse(env$cocoa_pct == "May-20",12.5,env$cocoa_pct)
env$cocoa_pct = as.numeric(env$cocoa_pct)
env$cocoa_pct = ifelse(env$cocoa_pct >1, env$cocoa_pct/100, env$cocoa_pct)

unique(env$grass_pct)
env$grass_pct = ifelse(env$grass_pct == "5-20",12.5,env$grass_pct)
env$grass_pct = ifelse(env$grass_pct == "May-20",12.5,env$grass_pct)
env$grass_pct = as.numeric(env$grass_pct)
env$grass_pct = ifelse(env$grass_pct >1, env$grass_pct/100, env$grass_pct)

unique(env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m == "10-20",15,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m == "Oct-20",15,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m == "20-30",25,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m == "30-40",35,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="5-10",7.5,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="05-Oct",7.5,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="2-5",3.5,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="02-May",3.5,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="1-2",1.5,env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m=="01-Feb",1.5,env$canopyHt_m)
env$canopyHt_m = as.numeric(env$canopyHt_m)
env$canopyHt_m = ifelse(env$canopyHt_m >1, env$canopyHt_m/100, env$canopyHt_m)

df[,c(2:17)] = env

####
### Add distance of trap to GRNP boundary : This was calculated externally in ArcGIS Pro
### Inside park = -ve distance
dist = read.csv("data2/Distance_to_GRNP_boundary.csv")
ord = match(x = df$SITE_ID, table = dist$IN_FID)
dist = dist[ord,]  
df$dist_f = dist

# write.csv(df, "data2/XY_No_Satellite_Gola.csv")



