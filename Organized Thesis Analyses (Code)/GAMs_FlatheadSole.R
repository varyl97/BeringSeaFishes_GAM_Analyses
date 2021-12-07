
# EGGS: Spawning Behavior Flathead Sole -----------------------------------
#the following code creates generalized additive models for eggs and larvae of flathead sole. 
#these analyses form the basis of my MS thesis. 
#egg data uses an averaged sea surface temperature for the month of March in the Southeastern Bering Sea
        #March index was chosen because it is two months before the peak of flathead egg CPUE, and thus March 
        #conditions are likely more relevant to spawning behavior than temperatures in later months. 

# Load in data ------------------------------------------------------------

fhsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhEggs.csv',header=TRUE,check.names=TRUE) #for egg GAMs
fhlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE) #for larval GAMs

##Load in regional temperature index for March (2 mos before peak egg CPUE in May), for egg GAMs 
reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(fhsub)){
  fhsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==fhsub$year[i]]}


# Spawning Behavior - Egg GAMs --------------------------------------------
#the following code generations GAMs that quantify variation in space and time of spawning behavior 
    #(as proxied by egg CPUE) in relation to regional temperature index close to typical time of spawning 

# Base Egg Model ----------------------------------------------------------

eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=fhsub,family=tw(link='log'),method='REML')

summary(eg.base)

windows(width=12,height=8)
plot(eg.base,shade=FALSE,page=1,
     seWithMean=TRUE,scale=0)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(eg.base,select=1,scheme=2,too.far=0.025,xlab='Longitude',
     ylab='Latitude',main='Base Egg Model, fh')
map("world",fill=T,col="snow4",add=T)
plot(eg.base,select=2,xlab='Day of Year',shade=FALSE,)
abline(h=0,col='sienna3',lty=2,lwd=2)


# Threshold Phenology Egg GAM ---------------------------------------------

temps<-sort(unique(reg.sst$SST)) #order by unique values of regional SST (March temps here)

bd<-4 #dictates essentially how many unique temperatures we check. can vary this depending on your system, how long you want the model to take, 
  #and probably some other factors that I don't know as well.. relatively arbitrary. smaller values for bd check more unique temperature values 

temps.in<-temps[bd:(length(temps)-bd)] #vector whose length is dependent on bd; this vector is what the for loop uses to test many models at different thresholds 

aic.pheno<-NA*(temps.in) #create an empty vector with the proper length to fill in later

thr.pheno<-as.list(1:(length(temps.in))) #another list with the same dimensions to fill in later

#this function below takes the unique temperature values and creates a model for each temperature as a threshold. 
#then, we use AIC to ask which model and at what threshold best explains variation in phenology? 
for(i in 1:length(temps.in)){
  fhsub$th<-factor(fhsub$reg.SST<=temps.in[i]) #how this is written is arbitrary, just sorts the reg.SST to below or equal to the given temp threshold
                                                #could also write it as above or equal to the given temp threshold 
                                                #as a result of this, for each row, the th column says either 'true' or 'false' depending on if the temperature
                                                          #is below or equal to the threshold temperature 
  thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+
                        s(lon,lat)+
                        s(bottom_depth,k=5)+
                        s(doy,by=th),
                      data=fhsub,family=tw(link='log'),method='REML') #standard gam formulation but the s(doy) term is where the threshold temp comes into play
  aic.pheno[i]<-AIC(thr.pheno[[i]]) #calculates AIC index for each model tested 
}

best.index.phe<-order(aic.pheno)[1] #now we're telling R to give us the model with the best (lowest) AIC value 
thr.pheno<-thr.pheno[[best.index.phe]]

windows()
plot(temps.in,aic.pheno,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.pheno)),
     main='Temperature Threshold Flexible Phenology',ylab='AIC Index',
     xlab='Temperature (degC)')
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.phe],lty=2,lwd=2,col='steelblue3') #this plot gives us a visual representation of every model tested and their respective AIC values 

summary(thr.pheno) #again, only care about the model with the best AIC value 

windows(width=12,height=8)
plot(thr.pheno,shade=FALSE,page=1,
     seWithMean=TRUE,scale=0)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(thr.pheno,select=4,main=paste('Below',round(temps.in[best.index.phe],digits=3)
                                                     ,sep=" "),
     shade=FALSE,seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(thr.pheno,select=3,main=paste('Above',round(temps.in[best.index.phe],digits=3),
                                                     sep=" "),
     shade=FALSE,seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)

windows()
par(mfrow=c(2,2))
gam.check(thr.pheno)


# Threshold Geography Egg GAM ---------------------------------------------

temps<-sort(unique(reg.sst$SST))
bd<-4 #change this to 4, more intermediate 
temps.in<-temps[bd:(length(temps)-bd)]

aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  fhsub$th<-factor(fhsub$reg.SST<=temps.in[i])
  thr.geo[[i]]<-gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                      s(lon,lat,by=th),data=fhsub,
                    family=tw(link='log'),method='REML')
  aic.geo[i]<-AIC(thr.geo[[i]])
} #add TH into grid.extent for true false based on condition 

best.index.geo<-order(aic.geo)[1]
thr.geo<-thr.geo[[best.index.geo]]
windows()
plot(temps.in,aic.geo,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.geo)),
     main='Temperature Threshold Flex Geography',xlab="Temperature (degC)")
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.geo],lty=2,lwd=2,col='steelblue3')

summary(thr.geo)

windows(width=12,height=8)
plot(thr.geo,page=1,scale=0,shade=FALSE,
     seWithMean=TRUE,scheme=2)

windows(width=12,height=8)
par(mfrow=c(1,2),oma=c(1,1,1,5))
plot(thr.geo,select=4,scheme=2,too.far=0.025,
     main=paste('Below',round(temps.in[best.index.geo],digits=3),'C',sep=" "),
     shade=FALSE,seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)
plot(thr.geo,select=3,scheme=2,too.far=0.025,
     main=paste('Above',round(temps.in[best.index.geo],digits=3),'C',sep=" "),
     shade=FALSE,seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)
gradientLegend(c(0,-5),color=c('red','orange','yellow'),ncol=5,side=4)

windows()
par(mfrow=c(2,2))
gam.check(thr.geo)


# Variable-Coefficient Phenology Egg GAM  ---------------------------------

vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                s(doy,by=reg.SST),data=fhsub,family=tw(link='log'),
              method='REML')
summary(vc.pheno)

windows(width=12,height=8)
plot(vc.pheno,shade=FALSE,
     page=1,scale=0,main='V-C Temp Flexible Phenology, fh Eggs',
     seWithMean=TRUE)

windows(width=16,height=8)
par(mfrow=c(1,2))
plot(vc.pheno,select=2,shade=FALSE,
     main='V-C Regional Temp, Flexible Phenology',xlab='Day of Year',
     seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(vc.pheno,select=4,shade=FALSE,
     main='V-C Reg. Temp, Deviation from Avg. Pheno Variation',
     xlab='Day of Year',seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)

windows()
par(mfrow=c(2,2))
gam.check(vc.pheno)


# Variable-Coefficient Geography Egg GAM ----------------------------------

vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
              s(lon,lat,by=reg.SST),data=fhsub,family=tw(link='log'),
            method='REML')
summary(vc.geo)

windows(width=12,height=8)
plot(vc.geo,shade=FALSE,
     page=1,seWithMean=TRUE,main='V-C Flex Geo')

windows(width=15,height=8)
par(mfrow=c(1,2))
plot(vc.geo,select=1,scheme=2,too.far=0.025,shade=FALSE,
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude',
     main='V-C fh Egg Flex Geo, Avg. Variation')
map("world",fill=T,col="snow4",add=T)
plot(vc.geo,select=4,scheme=2,too.far=0.025,shade=FALSE,
     xlab='Longitude',ylab='Latitude',seWithMean=TRUE,
     main='V-C Flex Geo, Deviation from Avg. Variation')
map("world",fill=T,col="snow4",add=T)

windows()
par(mfrow=c(2,2))
gam.check(vc.geo)


# Collation of all Egg Models  --------------------------------------------

eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=fhsub,family=tw(link='log'),method='REML')

summary(eg.base)

#threshold phenology
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]
aic.pheno<-NA*(temps.in)
thr.pheno<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  fhsub$th<-factor(fhsub$reg.SST<=temps.in[i])
  thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(bottom_depth,k=5)+s(doy,by=th),
                      data=fhsub,family=tw(link='log'),method='REML')
  aic.pheno[i]<-AIC(thr.pheno[[i]])
}

best.index.phe<-order(aic.pheno)[1]
thr.pheno<-thr.pheno[[best.index.phe]]
summary(thr.pheno)

#threshold geography
temps<-sort(unique(reg.sst$SST))
bd<-4 #change this to 4, more intermediate (checks more temperatures than bd = 10)
temps.in<-temps[bd:(length(temps)-bd)]

aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  fhsub$th<-factor(fhsub$reg.SST<=temps.in[i])
  thr.geo[[i]]<-gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                      s(lon,lat,by=th),data=fhsub,
                    family=tw(link='log'),method='REML')
  aic.geo[i]<-AIC(thr.geo[[i]])
} #add TH into grid.extent for true false based on condition 

best.index.geo<-order(aic.geo)[1]
thr.geo<-thr.geo[[best.index.geo]]
summary(thr.geo)

vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                s(doy,by=reg.SST),data=fhsub,family=tw(link='log'),
              method='REML')
summary(vc.pheno)

vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
              s(lon,lat,by=reg.SST),data=fhsub,family=tw(link='log'),
            method='REML')
summary(vc.geo)

# Saved models for loading later! -----------------------------------------
#MAKE SURE FOR THRESHOLDS that theyre either renamed as thr.geo or you're saving "thr.geo[[best.index.geo]]"
saveRDS(eg.base,file="./GAM Models/fh_egg_base.rds")
saveRDS(thr.pheno,file="./GAM Models/fh_egg_thr_pheno.rds")
saveRDS(temps.in,file="./GAM Models/fh_egg_temps_in_pheno.rds")
saveRDS(aic.pheno,file="./GAM Models/fh_egg_aic_pheno.rds")
saveRDS(best.index.phe,file="./GAM Models/fh_egg_best_index_pheno.rds")
saveRDS(thr.geo,file="./GAM Models/fh_egg_thr_geo.rds")
saveRDS(temps.in,file="./GAM Models/fh_egg_temps_in_geo.rds")
saveRDS(aic.geo,file="./GAM Models/fh_egg_aic_geo.rds")
saveRDS(best.index.geo,file="./GAM Models/fh_egg_best_index_geo.rds")
saveRDS(vc.pheno,file="./GAM Models/fh_egg_vc_pheno.rds")
saveRDS(vc.geo,file="./GAM Models/fh_egg_vc_geo.rds")

eg.base<-readRDS("./GAM Models/fh_egg_base.rds")
thr.pheno<-readRDS("./GAM Models/fh_egg_thr_pheno.rds")
temps.in<-readRDS("./GAM Models/fh_egg_temps_in_pheno.rds") 
      #these are based off of unique values in regional index, so "temps.in" doesn't change depending on model 
best.index.phe<-readRDS("./GAM Models/fh_egg_best_index_pheno.rds")
thr.geo<-readRDS("./GAM Models/fh_egg_thr_geo.rds")
best.index.geo<-readRDS("./GAM Models/fh_egg_best_index_geo.rds")
vc.pheno<-readRDS("./GAM Models/fh_egg_vc_pheno.rds")
vc.geo<-readRDS("./GAM Models/fh_egg_vc_geo.rds")


# Use Akaike Information Criterion to Check Model Performance -----------

aic.base<-AIC(eg.base)
aic.thrph<-AIC(thr.pheno)
aic.thrge<-AIC(thr.geo)
aic.vcph<-AIC(vc.pheno)
aic.vcgeo<-AIC(vc.geo)

aic.fhegg<-data.frame('model'=c('Base','Threshold Pheno','Threshold Geo',
                                'VC Pheno','VC Geo'),
                      'AIC_value'=c(aic.base,aic.thrph,aic.thrge,
                                    aic.vcph,aic.vcgeo))

windows()
plot(c(1:5),aic.fhegg$AIC_value,main='AIC Results for fh Egg Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.fhegg$AIC_value,labels=round(aic.fhegg$AIC_value),pos=c(4,3,3,3,2))
legend("topright",legend=c('Base','Threshold Pheno','Threshold Geo',
                             'VC Pheno','VC Geo'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)


# Larval Biogeography GAMs ------------------------------------------------
#the following code evaluates larval biogeography in relation to temperature and salinity
#temp and salinity measurements were taken in situ via conductivity-temperature-depth instruments 
      #on the same cruises where larval data were collected 

# Base Larval Biogeography Model ------------------------------------------

lv.base<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
               s(bottom_depth,k=5),
             data=fhlarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(lv.base,select=1,shade=FALSE,
     seWithMean=TRUE,scale=0,main='Base Larval Presence GAM, W/O Residuals')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.base,select=2,scheme=2,too.far=0.025,
     shade=FALSE,
     seWithMean=TRUE,scale=0)
map("world",fill=T,col="snow4",add=T)



# Larval Biogeography with Additive Salinity Interaction ------------------

lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                  s(bottom_depth,k=5)+
                  s(salinity),data=fhlarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

windows(width=12,height=8)
plot(lv.add.sal,page=1,scale=0,shade=FALSE,
     seWithMean=TRUE)

windows(width=14,height=8)
par(mfrow=c(1,2))
plot(lv.add.sal,select=1,seWithMean=TRUE,shade=FALSE,
     main='Seasonal Presence, Added Sal Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.add.sal,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows(width=12,height=8)
plot(lv.add.sal,select=4,shade=FALSE,
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Salinity',
     xlab='Salinity (PSU)')
abline(h=0,col='sienna3',lty=2,lwd=2)


# Larval Biogeography with Additive Temperature Interaction ---------------

lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature),data=fhlarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

windows()
plot(lv.add.temp,page=1,shade=FALSE,
     seWithMean=TRUE,main='Larval Log(Cper10m2+1) w Temp',scale=0)

windows(width=14,height=8)
par(mfrow=c(1,2))
plot(lv.add.temp,select=1,seWithMean=TRUE,shade=FALSE,
     main='Seasonal Presence, Added Temp Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.add.temp,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows(width=12,height=8)
plot(lv.add.temp,select=4,shade=FALSE,
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Temperature',
     xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)


# Biogeography with Individual Additive Temperature and Salinity Interactions --------

lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature)+s(salinity),data=fhlarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

windows()
plot(lv.temp.sal,page=1,shade=FALSE,
     main='Larval Log Presence Temp and Sal',
     seWithMean=TRUE,scale=0)

windows()
par(mfrow=c(2,2))
plot(lv.temp.sal,select=1,seWithMean=TRUE,shade=FALSE,
     main='Temp and Sal Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)
plot(lv.temp.sal,select=4,shade=FALSE,
     seWithMean=TRUE,main='Effect of Temp',xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=5,shade=FALSE,
     seWithMean=TRUE,main='Effect of Salinity',xlab='Salinity (psu)')
abline(h=0,col='sienna3',lty=2,lwd=2)


# Biogeography with Two-Dimensional Interaction between Temperature and Salinity --------

lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
             s(salinity,temperature),data=fhlarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)

windows()
plot(lv.2d,page=1,shade=FALSE,
     main='Larval Log Presence, 2D Temp and Sal',
     seWithMean=TRUE,scale=0,scheme=2)

windows()
par(mfrow=c(1,2))
plot(lv.2d,select=2,seWithMean=TRUE,shade=FALSE,
     main='Seasonal Presence, 2D Temp+Sal Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.2d,select=1,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
plot(lv.2d,select=4,scheme=2,main='Larval Log Presence, 2D Temp and Sal Effect',
     too.far=0.025,
     xlab='Temperature (degC)',ylab='Salinity (psu)')


# Collation of all Larval GAMs --------------------------------------------

lv.base<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
               s(bottom_depth,k=5),
             data=fhlarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                  s(bottom_depth,k=5)+
                  s(salinity),data=fhlarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature),data=fhlarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature)+s(salinity),data=fhlarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
             s(salinity,temperature),data=fhlarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)


# SAVE MODELS for reloading later -----------------------------------------

saveRDS(lv.base,file="./GAM Models/fh_larv_base.rds")
saveRDS(lv.add.sal,file="./GAM Models/fh_larv_addsal.rds")
saveRDS(lv.add.temp,file="./GAM Models/fh_larv_addtemp.rds")
saveRDS(lv.temp.sal,file="./GAM Models/fh_larv_tempsal.rds")
saveRDS(lv.2d,file="./GAM Models/fh_larv_2d.rds")

lv.base<-readRDS("./GAM Models/fh_larv_base.rds")
lv.add.sal<-readRDS("./GAM Models/fh_larv_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/fh_larv_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/fh_larv_tempsal.rds")
lv.2d<-readRDS("./GAM Models/fh_larv_2d.rds")


# Model Validation by Akaike Information Criterion  -----------------------

aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.tempsal<-AIC(lv.temp.sal)
aic.2d<-AIC(lv.2d)

aic.fhlarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Add Sal, Temp','Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.tempsal,aic.2d))

windows()
plot(c(1:5),aic.fhlarv$AIC_value,main='AIC Results for fh Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.fhlarv$AIC_value,labels=round(aic.fhlarv$AIC_value),pos=c(4,3,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Add Sal, Temp','Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)

base.test.lv<-logLik.gam(lv.base)
sal.test.lv<-logLik.gam(lv.add.sal)
temp.test.lv<-logLik.gam(lv.add.temp)
tempsal.test.lv<-logLik.gam(lv.temp.sal)
lv.test.2d<-logLik.gam(lv.2d)

test.scores<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                  'Add Sal, Temp','Sal-Temp 2D'),
                        'logLik_Value'=c(base.test.lv,sal.test.lv,
                                         temp.test.lv,tempsal.test.lv,
                                         lv.test.2d))
windows()
plot(c(1:5),test.scores$logLik_Value,main='Test LogLik Model Validation',
     col=c("#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='LogLik Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),test.scores$logLik_Value,labels=round(test.scores$logLik_Value),
     pos=c(4,3,3,3,2))
legend("bottomright",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Add Sal, Temp','Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)


# Salinity-Temp Hotspot Figures -------------------------------------------
#these figures represent where ideal temperature and salinity conditions are located in the SEBS 
    #and how catch of larvae interacts with those "hot spots" 
    #temp and salinity "hot spots" are determined by the best model output, the 2D temp-salinity model. 

viridis<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF", "#25858EFF", "#1E9B8AFF",
                            "#2BB07FFF", "#51C56AFF", "#85D54AFF", "#C2DF23FF", "#FDE725FF")) #viridis palette


stsweet<-fhlarv.ctd[fhlarv.ctd$temperature>1&fhlarv.ctd$temperature<10&
                      fhlarv.ctd$salinity<33&fhlarv.ctd$salinity>30,]#dim 249, 27
tt<-table(stsweet$year)

stsweet<-subset(stsweet,year%in%names(tt[tt>5]))
table(stsweet$year) #identify years with these qualities 
stsweet<-subset(stsweet,Cper10m2>0)

stsweet$s_color<-NA #s for stations
stsweet$month_nm<-NA
stsweet$c_color<-NA #c for catch 
stsweet<-stsweet %>% mutate(s_color=case_when(
  stsweet$month==4~'#482173FF',
  stsweet$month==5~'#2D708EFF',
  stsweet$month==6~'#2BB07FFF',
  stsweet$month==7~'#85D54AFF',
  stsweet$month==9~'#FDE725FF'))
stsweet<-stsweet %>%mutate(month_nm=case_when(
  stsweet$month==4~'April',
  stsweet$month==5~'May',
  stsweet$month==6~'June',
  stsweet$month==7~'July',
  stsweet$month==9~'September'))
stsweet<-stsweet%>%mutate(c_color=case_when(
  stsweet$month==4~'#482173FF',
  stsweet$month==5~'#2D708EFF',
  stsweet$month==6~'#2BB07FFF',
  stsweet$month==7~'#85D54AFF',
  stsweet$month==9~'#FDE725FF'
))

stsweet$c_color<-adjust_transparency(stsweet$c_color,alpha=0.2)#so bubbles don't cloud everything else
no.col<-adjust_transparency('hotpink3',alpha=0.3)

years<-sort(unique(stsweet$year))
tmp1<-1:ceiling(length(years)/4)
lg.text<-unique(stsweet$month_nm)
lg.values<-c('#482173FF','#2D708EFF','#2BB07FFF','#85D54AFF','#FDE725FF')
lg.order<-matrix(1:8,ncol=4,byrow=TRUE)

add_legend<-function(...){
  opar<-par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  on.exit(par(opar))
  plot(0,0,type='n',bty='n',xaxt='n',yaxt='n')
  legend(...)}

for(j in 1:length(tmp1)){ 
  windows(width=24,height=14)
  par(mfcol=c(2,2),omi=c(0.8,0.3,0.55,0.3),mai=c(0.2,0.4,0.4,0.1))
  for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
    plot(stsweet$lon[stsweet$year==years[i]],stsweet$lat[stsweet$year==years[i]],
         col=stsweet$s_color[stsweet$year==years[i]],pch=16,cex=2,
         main=as.character(years[i]),
         ylim=range(fhlarv.ctd$lat),xlim=range(fhlarv.ctd$lon),
         ylab=expression(paste("Latitude ("^0,'N)')),
         xlab=expression(paste("Longitude ("^0,'E)')))
    points(fhlarv.ctd$lon[fhlarv.ctd$Cper10m2==0&fhlarv.ctd$year==years[i]&
                            fhlarv.ctd$month==stsweet$month[i]],
           fhlarv.ctd$lat[fhlarv.ctd$Cper10m2==0&fhlarv.ctd$year==years[i]&
                            fhlarv.ctd$month==stsweet$month[i]],
           pch=4,cex=1.2,add=T,col=no.col)
    symbols(fhlarv.ctd$lon[fhlarv.ctd$Cper10m2>0&fhlarv.ctd$year==years[i]],
            fhlarv.ctd$lat[fhlarv.ctd$Cper10m2>0&fhlarv.ctd$year==years[i]],
            circles=log(fhlarv.ctd$Cper10m2[fhlarv.ctd$Cper10m2>0&fhlarv.ctd$year==years[i]]+1),
            inches=0.15,add=T,bg=(stsweet$c_color[stsweet$year==years[i]]))
    map("worldHires",fill=T,col="snow4",add=T)
  }
  add_legend("bottom",horiz=TRUE,legend=c('April','May','June','July','September','No Catch',
                                          'Proportional Log(Cper10m2+1)')[lg.order],
              col=c(lg.values,'hotpink3')[lg.order],lwd=c(3,3,3,3,3,NA,NA)[lg.order],lty=1,
              pch=c(NA,NA,NA,NA,NA,4,1)[lg.order],bg='lightgrey',bty='n',xpd=TRUE)
  
}


legend("bottomright",inset=c(-0.3,0),
       legend=c('April','May','June','July','September','No Catch'),
       col=c(lg.values,'hotpink3'),lwd=c(3,3,3,3,3,NA),lty=1,
       pch=c(NA,NA,NA,NA,NA,4),
       bg='lightgrey',xpd=T)
mtext("FH Sole Cper10m2 Catch Plotted Against T-S 'Sweet Spots'",outer=TRUE,cex=1,line=1)
