library(dplyr)
library(sjSDM)

# Load all Data sets
# 1 - Survey Data - Counts & Environment
survey = read.csv("data2/XY_No_Satellite_Gola.csv", row.names = 1)
rownames(survey) = 1:nrow(survey)

# Create Count Data
y = survey[,c(20:ncol(survey))]
y[y>0]=1
# As per Pichler 2021 - Remove OTU's with fewer than 3 occurrences
# He Removes Sites with fewer than 4 OTU's - But I am keeping all sites that have any insects
df_minOcc = survey
r = which(rowSums(y)>0) # Rows (trap sites with data) to keep
df_minOcc = df_minOcc[c(r),]
c = which(colSums(y)>2) # Columns to keep
c = c+19
df_minOcc = df_minOcc[,c(1:19,c)]
y = df_minOcc[,c(20:ncol(df_minOcc))]
y[y>0]=1

## Canonical Components
buffer = 30 # Using 30m as that provided highest explanatory AUC  

files = paste0("data2/sgdm/components/Can_Cpnts_k10_lsat_glcm-bandslsat_glcm-indices_", buffer,".csv")
canon = read.csv(files)
names(canon)[4:13] = c('cc1','cc2','cc3','cc4','cc5','cc6','cc7','cc8','cc9','cc10')

files = paste0("data2/sgdm/components/Can_Cpnts_k1_RADAR_", buffer,".csv")
radar = read.csv(files)
radar = radar[,4]

hansP = read.csv("data2/golaTRAPpoints_HansenYear_maxMin_30_500m.csv")
hans = hansP[,grep(buffer, colnames(hansP))]
hans = hans[c(r),] # Remove the rows that are not in the count data
# Set the Hansen data to the Max value of buffer level
hans=hans[,1]

distance = df_minOcc$distance_f

# Create DF that will be used to Train Model & Test Predictions
df_used = cbind(canon,hans,distance,radar, y)

# Spatial terms: for linear spatial model
XY = df_used[,c(2,3)]
names(XY) = c("X1", "X2")

### Run Model On Full Data Set
## Linear Model - Parameters suggested by Max Pichler
model.train = sjSDM(Y = as.matrix(y),
                    env = as.matrix(scale(df_used[,4:16]), lambda = 0.01, alpha = 0.2 ),
                    biotic = bioticStruct(lambda = 0.01, 
                                          alpha = 0.2, reg_on_Cov = FALSE, df = ncol(y)),
                    spatial = linear(scale(XY), ~0+X1+X2:X1:X2+I(X1^2)+I(X2^2), 
                                     lambda = 0.01, alpha = 0.2),
                    control = sjSDMControl(lr_reduce_factor = 0.9, scheduler = 10L), # improve convergence of optimization, don't change
                    learning_rate = 0.01, # don't change 
                    step_size = 30L, # don't change
                    family = binomial("probit"),
                    iter = 500, # Number of optimization steps, don't change
                    device = "gpu", 
                    sampling = 5000L, # Number of Monte-Carlo samples 
                    seed = 7)

# save(model.train, file = "data2/models/Final_K10.rad.hans.dist.Rdata")
###########################################################################

# load("data2/models/Final_K10.rad.hans.dist.Rdata")
ava = anova(model.train)
# save(ava, file = "data2/models/Final_K10.rad.hans.dist_ANOVA.Rdata")
# load("data2/models/Final_K10.rad.hans.dist_ANOVA.Rdata")
plot(ava, internal=TRUE)
ava$to_print
ava$results

# ava$species$R2_McFadden$Full
mean(ava$species$R2_McFadden_shared$R2) # Mean Species McF R2
mean(ava$sites$R2_McFadden_shared$R2) # Mean Site McF R2

part = plot(ava, internal = T)
sit = part$data$Sites
sum(part$data$Species$r2)
mean(part$data$Sites$r2)
mean(part$data$Species$r2)

### Get Variance Partition Values for Env, Space, Co-Distribution
vp = part$data$Species
vp$sp = rownames(vp)
vp[,1:3][vp[,1:3]>1] = 0
mean(vp$env)  
mean(vp$spa)    
mean(vp$codist) 

