# Model Validation Techniques ---------------------------------------------
#This script collates methods for testing the validity and robustness of GAM models.  

# Checking Variable Co-variation (Semi-Quantitative) ----------------------

#This section checks the co-variation of all sets of variables included in GAMs for all species. 
#Though not fully quantitative methods, I plan to utilize these methods and biological understanding of the system 
      #to validate the inclusion of variables within my models. 
#This code also is not complete. Not all study species are included (but will be eventually). 

library(corrplot)

# Loading Data ------------------------------------------------------------
apsub<-read.csv("./Ichthyo Data/Cleaned_Cut_ApEggs.csv")
fhsub<-read.csv("./Ichthyo Data/Cleaned_Cut_FhEggs.csv")
pksub<-read.csv("./Ichthyo Data/Cleaned_Cut_PkEggs.csv")
rxsub<-read.csv("./Ichthyo Data/Cleaned_Cut_RxEggs.csv")

aplarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv")
fhlarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv")
pklarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv")
yflarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv")
pclarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv")
nrslarv.ctd<-read.csv("./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv")

reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(apsub)){
  apsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==apsub$year[i]]} #these three species need March temperatures

for(i in 1:nrow(fhsub)){
  fhsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==fhsub$year[i]]}

for(i in 1:nrow(pksub)){
  pksub$reg.SST[i]<-reg.sst$SST[reg.sst$year==pksub$year[i]]}

for(i in 1:nrow(rxsub)){
  rxsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==rxsub$year[i]]}

reg.sst<-read.csv('./Environmental Data/May_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(yfsub)){
  yfsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==yfsub$year[i]]}


# Simple Data Investigation  ----------------------------------------------
#Compare raw Cper10m2 with log-transformed (log(Cper10m2+1)) data: 
#eggs: 
apsub$logcpue<-log(apsub$Cper10m2+1)
fhsub$logcpue<-log(fhsub$Cper10m2+1)
pksub$logcpue<-log(pksub$Cper10m2+1)
rxsub$logcpue<-log(rxsub$Cper10m2+1)

apraw<-ggplot(apsub,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Plaice")
fhraw<-ggplot(fhsub,aes(x=Cper10m2))+geom_histogram(color="black",fill='white',bins=40)+
  ggtitle("Flathead")
pkraw<-ggplot(pksub,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Pollock")
rxraw<-ggplot(rxsub,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Rex")
aptran<-ggplot(apsub,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Plaice, Transformed")+labs(x='log(Cper10m2+1)')
fhtran<-ggplot(fhsub,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Flathead, Transformed")+labs(x='log(Cper10m2+1)')
pktran<-ggplot(pksub,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Pollock, Transformed")+labs(x='log(Cper10m2+1)')
rxtran<-ggplot(rxsub,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Rex, Transformed")+labs(x='log(Cper10m2+1)')

invesall<-ggarrange(apraw,aptran,pkraw,pktran,fhraw,fhtran,rxraw,rxtran,ncol=2,nrow=4)+
  labs(title='Egg Data')+theme(plot.title=element_text(size=rel(1)))
windows(height=20,width=15)
invesall

aplarv.ctd$logcpue<-log(aplarv.ctd$Cper10m2+1)
fhlarv.ctd$logcpue<-log(fhlarv.ctd$Cper10m2+1)
pclarv.ctd$logcpue<-log(pclarv.ctd$Cper10m2+1)
pklarv.ctd$logcpue<-log(pklarv.ctd$Cper10m2+1)
yflarv.ctd$logcpue<-log(yflarv.ctd$Cper10m2+1)
nrslarv.ctd$logcpue<-log(nrslarv.ctd$Cper10m2+1)


apraw<-ggplot(aplarv.ctd,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Plaice")
fhraw<-ggplot(fhlarv.ctd,aes(x=Cper10m2))+geom_histogram(color="black",fill='white',bins=40)+
  ggtitle("Flathead")
pcraw<-ggplot(pclarv.ctd,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Cod")
pkraw<-ggplot(pklarv.ctd,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Pollock")
yfraw<-ggplot(yflarv.ctd,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Yellowfin")
nrsraw<-ggplot(nrslarv.ctd,aes(x=Cper10m2))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Northern Rock")
aptran<-ggplot(aplarv.ctd,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Plaice, Transformed")+labs(x='log(Cper10m2+1)')
fhtran<-ggplot(fhlarv.ctd,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Flathead, Transformed")+labs(x='log(Cper10m2+1)')
pctran<-ggplot(pclarv.ctd,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Cod, Transformed")
pktran<-ggplot(pksub,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Pollock, Transformed")+labs(x='log(Cper10m2+1)')
yftran<-ggplot(yflarv.ctd,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Yellowfin, Transformed")
nrstran<-ggplot(nrslarv.ctd,aes(x=logcpue))+geom_histogram(color='black',fill='white',bins=40)+
  ggtitle("Northern Rock, Transformed")

inveslarv<-ggarrange(apraw,aptran,fhraw,fhtran,pcraw,pctran,pkraw,pktran,yfraw,yftran,
                     nrsraw,nrstran,ncol=4,nrow=3)+labs(title='Larvae Data')+
  theme(plot.title=element_text(size=rel(1)))

windows(width=24,height=18)
inveslarv

# Test Co-Variation via Correlation Matrices --------------------------------------
#Egg gam variables: catch per 10m2, bottom depth, lat and lon, doy, reg.SST, year
#Trial 1: 
library(corrplot)

apcorr<-apsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
fhcorr<-fhsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
pkcorr<-pksub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]
rxcorr<-rxsub[c('bottom_depth','Cper10m2','year','doy','lat','lon','reg.SST')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pk<-cor(pkcorr)
rx<-cor(rxcorr)

windows(width=26,height=20)
par(mfrow=c(2,2),mar=c(1,1,1,0))
corrplot(ap,method="color",title='Plaice',tl.pos='tp')
corrplot(fh,method="color",title='Flathead',tl.pos='tp')
corrplot(pk,method="color",title='Pollock',tl.pos='tp')
corrplot(rx,method="color",title='Rex',tl.pos='tp')#bottom depth did not have values for any other variables 
                            #other correlations were pretty understandable (i.e. lon and lat negatively covaried)

#Trial 2: 
col<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF",
                        "#25858EFF", "#1E9B8AFF","#2BB07FFF", "#51C56AFF", "#85D54AFF", 
                        "#C2DF23FF", "#FDE725FF"))

apcorr<-apsub[c('logcpue','year','doy','reg.SST')]
fhcorr<-fhsub[c('logcpue','year','doy','reg.SST')]
pkcorr<-pksub[c('logcpue','year','doy','reg.SST')]
rxcorr<-rxsub[c('logcpue','year','doy','reg.SST')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pk<-cor(pkcorr)
rx<-cor(rxcorr)

windows(width=30,height=20)
par(mfrow=c(2,2),omi=c(0,0,1,0))
corrplot(ap,method="color",title="Alaska Plaice",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),addCoef.col='grey',cl.pos='n',type="upper")
corrplot(fh,method="color",title="Flathead Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),addCoef.col='grey',cl.pos='n',type="upper")
corrplot(pk,method="color",title="Pollock",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),addCoef.col='grey',cl.pos='n',type="upper")
corrplot(rx,method="color",title="Rex Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),addCoef.col='grey',cl.pos='n',type="upper")
mtext("Egg Data",outer=TRUE,line=3)

#Larvae: 
col<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF",
                        "#25858EFF", "#1E9B8AFF","#2BB07FFF", "#51C56AFF", "#85D54AFF", 
                        "#C2DF23FF", "#FDE725FF"))

apcorr<-aplarv.ctd[c('logcpue','year','doy','temperature','salinity')]
fhcorr<-fhlarv.ctd[c('logcpue','year','doy','temperature','salinity')]
pccorr<-pclarv.ctd[c('logcpue','year','doy','temperature','salinity')]
pkcorr<-pklarv.ctd[c('logcpue','year','doy','temperature','salinity')]
yfcorr<-yflarv.ctd[c('logcpue','year','doy','temperature','salinity')]
nrscorr<-nrslarv.ctd[c('logcpue','year','doy','temperature','salinity')]

ap<-cor(apcorr)
fh<-cor(fhcorr)
pc<-cor(pccorr)
pk<-cor(pkcorr)
yf<-cor(yfcorr)
nrs<-cor(nrscorr)

windows(width=28,height=20)
par(mfrow=c(2,3),omi=c(0,0,1,0))
corrplot(ap,method="color",title="Alaska Plaice",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
corrplot(fh,method="color",title="Flathead Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
corrplot(pc,method="color",title="Cod",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
corrplot(pk,method="color",title="Pollock",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
corrplot(yf,method="color",title="Yellowfin Sole",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
corrplot(nrs,method="color",title="Northern Rock",col=col(15),tl.col="black",
         mar=c(1,0,1,1.5),type="upper",addCoef.col='grey',cl.pos='n')
mtext('Larval Data',outer=TRUE,line=3)

# Compare Reduction in Mean Square Error  ---------------------------------
#This section looks at the reduction in MSE that occurs when models vary from the base formulation. 
#For eggs, variation means flexible phenology and/or geography and regional SST
#For larvae, variation means the inclusion of in situ temperature and salinity values 

#Flathead Sole:
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.179, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.097, smaller reduction than thr phenology 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.170

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.080 #geography models produce largest reduction in MSE

AIC(eg.base)-AIC(thr.pheno) #901.13

lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg #0.126

AIC(lv.base)-AIC(lv.2d) #381.28


#Alaska Plaice: 
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.035, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.149, larger reduction than thr phenology 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.022

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.135 #geography models produce largest reduction in MSE

var.second<-(summary(vc.geo)$scale-summary(thr.geo)$scale)/summary(thr.geo)$scale
var.second #0.017

AIC(eg.base)-AIC(thr.geo) #note change in AIC score - 463.75

lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg #0.151

lv.second<-(summary(lv.temp.sal)$scale-summary(lv.2d)$scale)/summary(lv.temp.sal)$scale
lv.second #0.070

AIC(lv.base)-AIC(lv.2d) #403.40

#Walleye Pollock: 
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.035, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.147, larger reduction than thr phenology 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.034

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.131 #geography models produce largest reduction in MSE

lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg #0.136

AIC(eg.base)-AIC(thr.geo) #620.57
AIC(lv.base)-AIC(lv.2d) #530.67

#Yellowfin Sole: 
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.111, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.086, larger reduction than thr phenology 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.042

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.032 #geography models produce largest reduction in MSE

lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg # + 0.194

AIC(eg.base)-AIC(thr.pheno) #507.03
AIC(lv.base)-AIC(lv.2d) #596.00

#Pacific Cod (larvae only): 
lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg #+ 0.094

aic.chg<-AIC(lv.base)-AIC(lv.2d)
aic.chg # +229.84

#Northern Rock Sole (larvae only): 
lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg #+ 0.118

aic.chg<-AIC(lv.base)-AIC(lv.2d)
aic.chg # + 344.99

#Rex Sole: (removed larval GAMs from thesis work due to lack of sampling/positive observations)
var.ratio.phe<-(summary(eg.base)$scale-summary(thr.pheno)$scale)/summary(eg.base)$scale
var.ratio.phe # positive difference of 0.033, meaning egg MSE was slightly larger than thr phenology model 

var.ratio.geo<-(summary(eg.base)$scale-summary(thr.geo)$scale)/summary(eg.base)$scale
var.ratio.geo # +0.188, largest reduction across all 

var.ratio.vcp<-(summary(eg.base)$scale-summary(vc.pheno)$scale)/summary(eg.base)$scale
var.ratio.vcp # +0.016

var.ratio.vcg<-(summary(eg.base)$scale-summary(vc.geo)$scale)/summary(eg.base)$scale
var.ratio.vcg # +0.101 #geography models produce largest reduction in MSE

lv.2d.chg<-(summary(lv.base)$scale-summary(lv.2d)$scale)/summary(lv.base)$scale
lv.2d.chg # + 0.052

thr.chg<-AIC(eg.base)-AIC(thr.geo)
thr.chg #647.63

lv.chg<-AIC(lv.base)-AIC(lv.2d)
lv.chg #113.18






























