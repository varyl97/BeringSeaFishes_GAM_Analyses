

# Model Validation Techniques ---------------------------------------------
#This script collates methods for testing the validity and robustness of GAM models.  

# Checking Variable Co-variation (Semi-Quantitative) ----------------------

#This section checks the co-variation of all sets of variables included in GAMs for all species. 
#Though not fully quantitative methods, I plan to utilize these methods and biological understanding of the system 
      #to validate the inclusion of variables within my models. 

library(corrplot)

# Loading Data ------------------------------------------------------------
apsub<-read.csv("./Ichthyo Data/Cleaned_Cut_ApEggs.csv")
fhsub<-read.csv("./Ichthyo Data/Cleaned_Cut_FhEggs.csv")
pksub<-read.csv("./Ichthyo Data/Cleaned_Cut_PkEggs.csv")
yfsub<-read.csv("./Ichthyo Data/Cleaned_Cut_YfEggs.csv")

aplarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv")
fhlarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv")
pklarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv")
yflarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv")

reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(apsub)){
  apsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==apsub$year[i]]} #these three species need March temperatures

for(i in 1:nrow(fhsub)){
  fhsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==fhsub$year[i]]}

for(i in 1:nrow(pksub)){
  pksub$reg.SST[i]<-reg.sst$SST[reg.sst$year==pksub$year[i]]}

reg.sst<-read.csv('./Environmental Data/May_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(yfsub)){
  yfsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==yfsub$year[i]]}

# Test Co-Variation via Scatter Plots --------------------------------------
#Egg gam variables: catch per 10m2, bottom depth, lat and lon, doy, reg.SST, year
#Trial 1: 
apcorr<-apsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
fhcorr<-fhsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
pkcorr<-pksub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
yfcorr<-yfsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pk<-cor(pkcorr)
yf<-cor(yfcorr)

windows(width=20,height=20)
par(mfrow=c(2,2))
corrplot(ap,method="color")
corrplot(fh,method="color")
corrplot(pk,method="color")
corrplot(yf,method="color") #bottom depth did not have values for any other variables 
                            #other correlations were pretty understandable (i.e. lon and lat negatively covaried)

#Trial 2: 
col<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF",
                        "#25858EFF", "#1E9B8AFF","#2BB07FFF", "#51C56AFF", "#85D54AFF", 
                        "#C2DF23FF", "#FDE725FF"))

apcorr<-apsub[c('Cper10m2','year','doy','reg.SST')]
fhcorr<-fhsub[c('Cper10m2','year','doy','reg.SST')]
pkcorr<-pksub[c('Cper10m2','year','doy','reg.SST')]
yfcorr<-yfsub[c('Cper10m2','year','doy','reg.SST')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pk<-cor(pkcorr)
yf<-cor(yfcorr)

windows(width=20,height=20)
par(mfrow=c(2,2))
corrplot(ap,method="color",title="Alaska Plaice",col=col(15),tl.col="black",
         mar=c(0,0,1,0))
corrplot(fh,method="color",title="Flathead Sole",col=col(15),tl.col="black",
         mar=c(0,0,1,0))
corrplot(pk,method="color",title="Pollock",col=col(15),tl.col="black",
         mar=c(0,0,1,0))
corrplot(yf,method="color",title="Yellowfin Sole",col=col(15),tl.col="black",
         mar=c(0,0,1,0))

#Larvae: 
col<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF",
                        "#25858EFF", "#1E9B8AFF","#2BB07FFF", "#51C56AFF", "#85D54AFF", 
                        "#C2DF23FF", "#FDE725FF"))

apcorr<-aplarv.ctd[c('Cper10m2','year','doy','temperature','salinity')]
fhcorr<-fhlarv.ctd[c('Cper10m2','year','doy','temperature','salinity')]
pkcorr<-pklarv.ctd[c('Cper10m2','year','doy','temperature','salinity')]
yfcorr<-yflarv.ctd[c('Cper10m2','year','doy','temperature','salinity')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pk<-cor(pkcorr)
yf<-cor(yfcorr)

windows(width=20,height=20)
par(mfrow=c(2,2))
corrplot(ap,method="color",title="Alaska Plaice",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5))
corrplot(fh,method="color",title="Flathead Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5))
corrplot(pk,method="color",title="Pollock",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5))
corrplot(yf,method="color",title="Yellowfin Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5))



# Compare Reduction in Mean Square Error  ---------------------------------
#This section looks at the reduction in MSE that occurs when models vary from the base formulation. 
#For eggs, variation means flexible phenology and/or geography and regional SST
#For larvae, variation means the inclusion of in situ temperature and salinity values 

#Flathead Sole:
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.06, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.167, larger reduction than thr phenology 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.045

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.136 #geography models produce largest reduction in MSE

#Yellowfin Sole: 











































