library(dplyr)
library(sjSDM)
library(glue)

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

## Canonical Components
buffer = 30 # Using 30m as trials showed it provided highest explanatory AUC  

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

owk = list() # To store train & test AUC values from 5 model runs using best parameters  

for (k in 1:4) {  #### Run model 4 times (each one splits data into 5 folds for x-validation) to get mean AUC
# k=1  
ac.mn = c()  # Individual model Predictive (test) AUC
ac.exp = c() # Individual model Explanatory AUC
sp = c()     # Species Names per AUC
  
#####################################
# Split data into Training & Testing
df_used$fold <- sample(factor(rep(1:5, length.out=nrow(df_used)), labels=paste0("Fold", 1:5))) # Add number to split data into 5 folds
  
for (kk in 1:5) { # For each one of the 5 folds
# kk=2 
  
train  <- df_used[df_used$fold != paste0("Fold", kk), -ncol(df_used)] # Split data into Training
test   <- df_used[df_used$fold == paste0("Fold", kk), -ncol(df_used)] # and test
    
# Create Count Data
y.train = train[,c(17:ncol(train))]
y.train[y.train>0]=1
# As per Pichler 2021 - Remove OTU's with fewer than 3 occurrences
# He Removes Sites with fewer than 4 OTU's - But I am keeping all sites that have any insects
train.minOcc = train
r = which(rowSums(y.train)>0) # Rows to keep
train.minOcc = train.minOcc[c(r),]
c = which(colSums(y.train)>2) # Columns to keep
train.minOcc = train.minOcc[,c(1:16,c+16)]
y.train = train.minOcc[,c(17:ncol(train.minOcc))]
    
y.test = test[,c(17:ncol(test))]
y.test = y.test[ ,colnames(y.test)%in%names(c)]
    
    
# Create Training and Testing Data
env.train = train.minOcc[,c(4:16)]
env.test = test[,c(4:16)]
    
# Spatial terms: for linear spatial model
XY.train = train.minOcc[,c(2:3)]
names(XY.train) = c("X1", "X2")
    
XY.test = test[,c(2,3)]
names(XY.test) = c("X1", "X2")
    
s.otu.train = as.matrix(y.train) %>% unname
s.otu.test = as.matrix(y.test) %>% unname
    
## Linear Model - Parameters suggested by Max Pichler
model.train = sjSDM(Y = as.matrix(y.train),
                    env = as.matrix(scale(env.train), lambda = 0.01, alpha = 0.2 ),
                    biotic = bioticStruct(lambda = 0.01, 
                                          alpha = 0.2, reg_on_Cov = FALSE, df = ncol(y)),
                    spatial = linear(scale(XY.train), ~0+X1+X2:X1:X2+I(X1^2)+I(X2^2), 
                                         lambda = 0.01, alpha = 0.2),
                    control = sjSDMControl(lr_reduce_factor = 0.9, scheduler = 10L), # improve convergence of optimization, don't change
                    learning_rate = 0.01, # don't change 
                    step_size = 30L, # don't change
                    family = binomial("probit"),
                    iter = 500, # Number of optimization steps, don't change
                    device = "gpu", 
                    sampling = 5000L, # Number of Monte-Carlo samples 
                    seed = 7)  
    
#############################
    ## add prediction(auc) data
auc.all = data.frame(otu = as.character(names(y.train)), auc.pred = 0.1, auc.exp = 0.1 )
  
for (pred in c('explain', 'pred')) {
# explanatory AUC
    if (pred == 'explain') {
        newdd = NULL
        newsp = NULL
        sp.df = y.train
        set = 'explain'
      }
      # predictive AUC
      if (pred == 'pred') {
        newdd = scale(env.test)
        newsp = scale(XY.test)
        sp.df = y.test
        set = 'test'
      }
      
model1 = model.train
sp.df1 = sp.df
      
pred.df = predict(model1, newdata = newdd, SP = newsp, type = 'link')
      
sp.df1 = data.frame(sp.df1)
      
## Add a row at the bottom which says if the species is counted of not in the given training split
sp.df1 = rbind(sp.df1, count = (base::colSums(sp.df1)>0 )*1 )
      
# some OTUs don't occur in test set
# Only keep predictions for species that occur in the test split
pred.df = pred.df[ , which(sp.df1[nrow(sp.df1),] == 1)]
sp.df1 = sp.df1[1:nrow(pred.df), which(sp.df1[nrow(sp.df1), ] == 1)]
      
# .... calculate AUC
met.df = sapply(1:ncol(sp.df1), function(i) Metrics::auc(sp.df1[,i], pred.df[,i]))
auc.mean = mean(met.df)
      
# ... make long table
if (pred=='explain') {
    auc.te = data.frame(auc.exp = met.df, otu = names(sp.df1) )
    } else if (pred=='pred') {
        auc.te = data.frame(auc.pred = met.df, otu = names(sp.df1) ) }
auc.all = left_join(auc.all, auc.te, by = 'otu', suffix = c('', glue('.{1}')), copy=T)
} 
    
auc.all = select(auc.all, -'auc.pred',-'auc.exp')
auc.all$incidence = colSums(y.train)
auc.all = auc.all %>% rename(auc.pred = 3, auc.train = 2) 
auc.all = auc.all[order(auc.all$otu),]
    
ac.mn = c(ac.mn, auc.all$auc.pred)
ac.exp = c(ac.exp, auc.all$auc.train)
sp = c(sp, auc.all$otu)
} # End of fold of KK folds  
##########################################################################
### to get a mean predictability per species per run of K
auc.mean = as.data.frame(cbind(sp, ac.mn))
auc.mean$ac.mn = as.numeric(auc.mean$ac.mn)
  
auc.exp = as.data.frame(cbind(sp, ac.exp))
auc.exp$ac.exp = as.numeric(auc.exp$ac.exp)
  
mm = data.frame(sp = unique(auc.mean$sp))
mx = data.frame(sp = unique(auc.exp$sp))
  
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
  
for (i in 1:length(mx$sp)) {
    # i=1
    mmx = which(mx$sp[i] == auc.exp$sp)
    for (ii in 1:length(mmx)) {
      # ii=1
      mx[i,ii+1] = auc.exp[mmx[ii],"ac.exp"]
    }
  }
mmmx = apply(X = mx[,2:ncol(mx)], MARGIN = 1, FUN = mean, na.rm=T)
mx$mean.auc = mmmx
mm$meanEXPauc = mmmx
  
mm$incidence = colSums(y)[match(mm$sp, colnames(y))]
mm = mm[order(mm$sp),]
  
table(is.nan(mm$mean.auc))
table(mm$mean.auc>=0.7)
table(mm$incidence)

## Add the current XVal data to a list
owk = append(owk, list(mm))
## Get the species names that have had a High predictive AUC in the current run & assign them to a variable 
assign(paste0("HiAUC",k), mm$sp[mm$mean.auc>=0.7 & !is.nan(mm$mean.auc)])
} ## End of Run of K runs
#### owk is a list of the 4 individual sets of 5-fold cross validation results
# save(owk, file = "../data2/Results_AUC_4Runs_5F-XVAL.Rdata")

## Combine all species names that have had a predictive AUC >=0.7
hiauc = c(HiAUC1, HiAUC2, HiAUC3, HiAUC4)#, HiAUC5, HiAUC6, HiAUC7)#, HiAUC8, HiAUC9) #, HiAUC10)
## How many species have a predictive AUC >=0.7 at least once in the 4 runs of Cross Validation??
length(unique(hiauc))
a = as.matrix(table(hiauc))
## How many species have a Predictive AUC>=0.7 in at least 2 of the 4 Cross-Validation runs??
length(which(a[,1]>=2))
## Get the names of those species that had at least 2 AUC of >=0.7
a = a[a[,1]>=2,]
n = names(a)

## Get Pred & Exp AUC values for each run of cross validation
for (j in 1:length(owk)) {
o = owk[[j]]
ex = o[,c(1,8)]
pd = o[,c(1,7)]
if (j==1) {ooPred = pd; ooExp = ex}
if (j!=1) {ooPred = full_join(ooPred,pd, by="sp"); ooExp = full_join(ooExp,ex, by="sp")}
}

colnames(ooPred)[2:5] = c("Pred_Run_1", "Pred_Run_2", "Pred_Run_3", "Pred_Run_4")
colnames(ooExp)[2:5] = c("Exp_Run_1", "Exp_Run_2", "Exp_Run_3", "Exp_Run_4")

## How many sp Had a predictive AUC >= 0.7 in at least 1 of the runs??
length(which(apply(ooPred[,2:5], 1, max)>=0.7))

## How many sp Had a mean predictive AUC >= 0.7??
length(which(apply(ooPred[,2:5], 1, mean)>=0.7))

## How many sp Had a mean predictive AUC >= 0.7 in at least 2 of the 4 (half of the) XVal runs??
oop = ooPred
oop[,2:5][oop[,2:5]>=0.7] = 1
oop[,2:5][oop[,2:5]<1] = 0
length(which(apply(oop[,2:5], MARGIN = 1, sum)>=2))
## Get species names that have a high predictive AUC in at least 2 of the XVal runs
n = oop$sp[apply(oop[,2:5], MARGIN = 1, sum)>=2]
n = n[!is.na(n)]
# length(n)
times.hi.auc = apply(oop[,2:5], MARGIN = 1, sum)
which(apply(oop[,2:5], MARGIN = 1, sum)>=2)

Exp_MEAN = apply(ooExp[,2:5], MARGIN = 1, FUN = mean, na.rm=T)
Pred_MEAN = apply(ooPred[,2:5], MARGIN = 1, FUN = mean, na.rm=T)

## Combine data to save
AUC = as.data.frame(cbind(ooExp,ooPred[,-1],Exp_MEAN,Pred_MEAN,times.hi.auc, incidence = o[,9]))
# write.csv(AUC, "data2/Results_AUC_perSp_Final_Model.csv",row.names = F)

## Get names of the species that had Hi AUC in at least 2 of the 4 runs
ac = AUC[AUC$times.hi.auc>=2 & !is.na(AUC$times.hi.auc),]
sp = ac$sp
# write.csv(sp, "data2/Results_Species_HiAuc_min2of4.csv", row.names = F)

#########################
### Test Correlation of AUC vs Prevalence
# AUC = read.csv("data2/Results_AUC_perSp_MaxesModels.csv")
cor.test(AUC$incidence, AUC$Pred_MEAN, use = "complete.obs", method = "pearson") # correlation between prevalence and Predictive AUC
cor.test(AUC$incidence, AUC$Exp_MEAN, use = "complete.obs", method = "pearson") # correlation between prevalence & Explanatory AUC
