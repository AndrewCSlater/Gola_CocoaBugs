library(conflicted)
library(tidyverse)
library(sjSDM)
library(glue)
library(pROC)


conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer('colSums', 'base')


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
r = which(rowSums(y)>0) # Rows to keep
df_minOcc = df_minOcc[c(r),]
c = which(colSums(y)>2) # Columns to keep
c = c+19
df_minOcc = df_minOcc[,c(1:19,c)]
y = df_minOcc[,c(20:ncol(df_minOcc))]
y[y>0]=1

## Create Environmental data
env_field = df_minOcc[,c(4, 7:17)] # Field measured env data, Land use & activity
# NA's are not permitted in any data
env_field[is.na(env_field)] = 0

# Spatial terms: for linear spatial model
XY = df_minOcc[,c(6,5)]
names(XY) = c("X1", "X2")

# Create DF that will be used to Train Model & Test Predictions
df_used = cbind(df_minOcc$SITE_ID, XY,env_field,y)
names(df_used)[1] = "Plot_ID"

owk = list() # To store train & test AUC values from 5 model runs using best parameters  
ac.mn = c()
sp = c()
for (kk in 1:5) { # Running model 5 times to get average performance
# kk=1  
  
#####################################
# Split data into Training & Testing
sample <- caTools::sample.split(df_used[,1], SplitRatio = 0.75)
train  <- subset(df_used, sample == TRUE)
test   <- subset(df_used, sample == FALSE)
  
## Create folds for testing
#Randomly shuffle the data
train<-train[sample(nrow(train)),]
#Create 5 equally size folds
fold.id <- cut(seq(1,nrow(train)),breaks=5,labels=FALSE)
train = cbind(fold.id, train)
train = train[order(train$Plot_ID),]

# Create Count Data
y.train = train[,c(16:ncol(train))]
y.train[y.train>0]=1
# As per Pichler 2021 - Remove OTU's with fewer than 3 occurrences
# He Removes Sites with fewer than 4 OTU's - But I am keeping all sites that have any insects
train.minOcc = train
r = which(rowSums(y.train)>0) # Rows to keep
train.minOcc = train.minOcc[c(r),]
c = which(colSums(y.train)>2) # Columns to keep
train.minOcc = train.minOcc[,c(1:15,c+15)]
y.train = train.minOcc[,c(1,16:ncol(train.minOcc))]
  
y.test = test[,c(15:ncol(test))]
y.test = y.test[ ,colnames(y.test)%in%names(c)]
  
  
# Create Training and Testing Data
env.train = train.minOcc[,c(1,6:15)] ### Change to 5:15 to Include distance to Park
env.test = test[,c(5:14)]
  
# Spatial terms: for linear spatial model
XY.train = train.minOcc[,c(1,3:4)]
names(XY.train) = c("fold.id", "X1", "X2")
  
XY.test = test[,c(2,3)]
names(XY.test) = c("X1", "X2")

s.otu.train = as.matrix(y.train[,-c(1:3)]) %>% unname
s.otu.test = as.matrix(y.test[,-c(1:2)]) %>% unname

i=1 # To get the parameter values out of a list of 1 


## Run a Linear Model
model.train = sjSDM(Y = as.matrix(s.otu.train),
                    env = as.matrix(scale(env.train[,-1]),
                    spatial = linear(scale(XY.train[,-1]), ~0+X1+X2:X1:X2+I(X1^2)+I(X2^2)),
                    iter = 200, device = "gpu")) ### NB: Use device = "cpu" if you can't run on "gpu" (an Nvidia graphics card) 


#############################
## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(y.train)[-c(1:3)]), auc.test = 0.1, auc.exp = 0.1 )

for (pred in c('explain', 'test')) {
  # explanatory AUC
  if (pred == 'explain') {
    newdd = NULL
    newsp = NULL
    otudd = s.otu.train
    set = 'explain'
  }
  # predictive AUC
  if (pred == 'test') {
    newdd = scale(env.test)
    newdd[is.nan(newdd)] = 0
    newsp = scale(XY.test)
    otudd = s.otu.test
    set = 'test'
  }
  
  model1 = model.train
  otudd1 = otudd
  
pred.dd = apply(abind::abind(lapply(1:3, function(i) predict(model1, newdata = newdd, SP = newsp, type = 'link')) , along = -1L), 2:3, mean) %>% unname
  
  otudd1 = data.frame(otudd1)
  names(otudd1) = names(y.train)[-c(1:3)]
  
  otudd1 = rbind(otudd1, count = (base::colSums(otudd1)>0 )*1 )
  
  # some OTUs don't occur in test set
  pred.dd = pred.dd[ , which(otudd1[nrow(otudd1),] == 1)]
  
  otudd1 = otudd1[1:nrow(pred.dd), which(otudd1[nrow(otudd1), ] == 1)]
  
  otudd.pa = (otudd1>0)*1
  
  # .... calculate AUC
  roc.dd = sapply(1:ncol(otudd1), function(j) as.numeric(pROC::roc( response = otudd.pa[,j], predictor = pred.dd[,j], direction = '<', quiet = T)$auc))
  
  auc.mean = mean(roc.dd)
  
  # ... make long table
  if (pred=='explain') {
    auc.te = data.frame(auc.exp = roc.dd, otu = names(otudd1) )
  } else if (pred=='test') {
    auc.te = data.frame(auc.test = roc.dd, otu = names(otudd1) ) }
  auc.all = left_join(auc.all, auc.te, by = 'otu', suffix = c('', glue('.{i}')), copy=T)
} 

auc.all = select(auc.all, -'auc.test',-'auc.exp')
auc.all$incidence = colSums(s.otu.train)
auc.all = auc.all %>% rename(auc.test = 3, auc.train = 2) 
auc.all = auc.all[order(auc.all$otu),]

ac.mn = c(ac.mn, auc.all$auc.test)
sp = c(sp, auc.all$otu)
owk = append(list(auc.all),owk)
}

##########################################################################
### This section used in run of 5 random data splits for training/testing
### to get a mean predictability per species
auc.mean = as.data.frame(cbind(sp, ac.mn))
auc.mean$ac.mn = as.numeric(auc.mean$ac.mn)

mm = data.frame(sp = unique(auc.mean$sp))

for (i in 1:length(mm$sp)) {
  # i=1
  mmm = which(mm$sp[i] == auc.mean$sp)
  for (ii in 1:length(mmm)) {
    # ii=1
    mm[i,ii+1] = auc.mean[mmm[ii],"ac.mn"]
  }
}
mmmm = apply(X = mm[,2:ncol(mm)], MARGIN = 1, FUN = mean, na.rm=T)
mm$mean.auc = mmmm

table(is.nan(mm$mean.auc))
table(mm$mean.auc>=0.7)

mm$incidence = colSums(y)[match(mm$sp, colnames(y))]

# write.csv(mm, "data2/Results_AUC_perSp_FieldData_model.csv",row.names = F)

auc_train = mean(auc.all$auc.train, na.rm = T) #
auc_mean = mean(auc.all$auc.test, na.rm = T)   # 

