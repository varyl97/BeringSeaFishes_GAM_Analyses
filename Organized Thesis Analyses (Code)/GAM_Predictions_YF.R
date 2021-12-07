##GAM Predictions: Yellowfin Sole
#the following code produces predictions and predictive figures based off the best model 
#for Yellowfin Sole, the best model is the threshold geography model wherein spawning geography 
#varies differently above and below 4.143 *C. The best phenology model was also the threshold formulation
#where the threshold temperature was 4.175 *C. 

# Preliminary Data Loading and Formatting ---------------------------------

#link yfsub to regional indices:
yfsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfEggs.csv',header=TRUE,check.names=TRUE) #for egg GAMs
yflarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv',header=TRUE,check.names=TRUE) #for larval GAMs

reg.sst<-read.csv('./Environmental Data/May_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(yfsub)){
  yfsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==yfsub$year[i]]}

#load GAMs
eg.base<-readRDS("./GAM Models/yf_egg_base.rds")
thr.pheno<-readRDS("./GAM Models/yf_egg_thr_pheno.rds")
best.index.phe<-readRDS("./GAM Models/yf_egg_best_index_phe.rds")
thr.geo<-readRDS("./GAM Models/yf_egg_thr_geo.rds")
vc.pheno<-readRDS("./GAM Models/yf_egg_vc_pheno.rds")
vc.geo<-readRDS("./GAM Models/yf_egg_vc_geo.rds") #only using eg.base, thr.pheno, and thr.geo below. 

lv.base<-readRDS("./GAM Models/yf_larvae_base.rds")
lv.add.sal<-readRDS("./GAM Models/yf_larvae_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/yf_larvae_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/yf_larvae_addtempsal.rds")
lv.2d<-readRDS("./GAM Models/yf_larvae_2d.rds") #only using lv.base and lv.2d in this code. 

#get map 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) 

# Visualize on a regular grid ---------------------------------------------
windowsFonts(A="Times New Roman") #for axes labels and plot titles, eventually

#visualize results of best model by predicting on a regular spaced grid: 
#prediction grid
nlat=120
nlon=120
latd=seq(min(yfsub$lat),max(yfsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(yfsub$lon),max(yfsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          yfsub$lat,yfsub$lon)
  grid.extent$dist[k]<-min(dist)
}

# Plot Base Geography Model -----------------------------------------------
grid.extent$year<-as.numeric(2008) #for base, just pick a year that has coverage; change this depending on species
grid.extent$doy<-as.numeric(median(yfsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(yfsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(yfsub$reg.SST,na.rm=TRUE) 
grid.extent$pred<-predict(eg.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(yfsub$lon),ylim=range(yfsub$lat),main='Yellowfin Sole Distribution, Eggs',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(yfsub$lon[yfsub$Cper10m2==0],yfsub$lat[yfsub$Cper10m2==0],pch='+',col='white')
symbols(yfsub$lon[yfsub$Cper10m2>0],
        yfsub$lat[yfsub$Cper10m2>0],
        circles=log(yfsub$Cper10m2+1)[yfsub$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

# Plot Base Phenology Model -----------------------------------------------
grid.extent2<-data.frame('lon'=rep(-155,100),
                         'lat'=rep(51,100),'doy'=seq(min(yfsub$doy),max(yfsub$doy),length=100),
                         'year'=rep(2008,100),'bottom_depth'=rep(median(yfsub$bottom_depth,na.rm=TRUE),100))#setting up another clean grid on which to predict
grid.extent2$pred<-predict(eg.base,newdata=grid.extent2)
grid.extent2$se<-predict(eg.base,newdata=grid.extent2,se=T)[[2]] #select the standard error 
grid.extent2$pred.u<-grid.extent2$pred+1.96*grid.extent2$se #calculate upper confidence interval (95% CI)
grid.extent2$pred.l<-grid.extent2$pred-1.96*grid.extent2$se #calculate lower CI 


windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
plot(grid.extent2$doy,grid.extent2$pred,main='Yellowfin Sole Base Phenology, Eggs',type='l',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),xlab='Day of Year',cex.lab=1,
     cex.axis=1,cex.main=1,cex.axis=0.9,xlim=range(yfsub$doy),
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred.l)),
     col='honeydew2',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col='honeydew2',lty=0)
lines(grid.extent2$doy,grid.extent2$pred,col='grey43')
legend('topleft',legend=c(expression(paste("(log(C/(10m"^2,')+1)')),'95% CI'),
       col=c('grey43','honeydew2'),pch=c(NA,15),lty=c(1,NA),lwd=2,cex=1)
abline(h=0,col='grey79',lty=2,lwd=1.5)

#with base graphics: 
windows()
par(mai=c(1,1,0.5,0.5))
plot(eg.base,select=2,main='Yellowfin Sole Base Phenology, Eggs',
     seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies (edf: 8.796)',ylim=c(-1,1))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)

#TEMP EFFECT: Calculate Differences Due to Different Temperature Regimes Based on Best Model --------
#start with threshold phenology model to find differences between two predictions (base prediction and this prediction) to calculate local slopes 
nlat=120
nlon=120
latd=seq(min(yfsub$lat),max(yfsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(yfsub$lon),max(yfsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          yfsub$lat,yfsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-2002
grid.extent$doy<-median(yfsub$doy)
grid.extent$reg.SST<-mean(yfsub$reg.SST) #threshold temp chosen by AIC values
grid.extent$th<-"TRUE"
grid.extent$bottom_depth<-median(yfsub$bottom_depth,na.rm=T)
grid.extent$pred<-predict(eg.base,newdata=grid.extent)
grid.extent$se<-predict(eg.base,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se #95% CI here
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA #remove predictions that are too far from positive data values
grid.extent$reg.SST<-mean(yfsub$reg.SST)
grid.extent$th<-"FALSE"
grid.extent$pred2<-predict(thr.pheno,newdata=grid.extent)
grid.extent$se2<-predict(thr.pheno,newdata=grid.extent,se=T)[[2]]
grid.extent$pred2.u<-grid.extent$pred2+1.96*grid.extent$se
grid.extent$pred2.l<-grid.extent$pred2-1.96*grid.extent$se
grid.extent$diff<-grid.extent$pred2-grid.extent$pred #calculate difference between two regimes

grid.extent$sig.pos<-c(grid.extent$pred2.l>grid.extent$pred.u) #isolate areas where there is a higher predicted CPUE at a higher temperature
grid.extent$sig.neg<-c(grid.extent$pred2.u<grid.extent$pred.l)
grid.extent$pos.diff<-grid.extent$diff*grid.extent$sig.pos #calculate areas with a significant positive difference at a higher temperature
grid.extent$neg.diff<-grid.extent$diff*grid.extent$sig.neg
max.slope<-max(grid.extent$diff,na.rm=T)

windows(width=15,height=15)
par(mai=c(1,1,0.5,0.5))
image.plot(lond,latd,t(matrix(grid.extent$diff,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=hcl.colors(100,"PRGn"),ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')), #PRGn diverges more clearly, helping interpretation
           xlim=range(yfsub$lon),ylim=range(yfsub$lat),main='Change in YFS Egg Distribution From Base to Threshold Pheno Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),
           legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)#would prefer to have legend within plot margins, and for all font to be times, but not sure how to do that. 
map("worldHires",fill=T,col="seashell2",add=T)

#Plot the two threshold phenology smooths together, one from the base model and one from the threshold geography model to see the temp effect: 
col<-adjustcolor('tomato4',alpha.f=0.3)

windows()
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5))
plot(thr.pheno,select=4,main='Yellowfin Sole Phenology, Eggs',seWithMean=TRUE,
     ylim=c(-3,3))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5),new=TRUE)
plot(thr.pheno,select=3,seWithMean=TRUE,shade=TRUE,shade.col=col,ylim=c(-3,3))
legend('topright',legend=c('Below','Above'),col=c(NA,col),lwd=c(2,2),cex=0.8)
mtext(c("Day of Year","Anomalies in log(CPUE+1)"),side=c(1,2),line=2.5)

#For added information: threshold phenology prediction (this was the best model to explain phenology):
grid.extent3<-data.frame('lon'=rep(-155,100),'lat'=rep(51,100),'doy'=seq(min(yfsub$doy),max(yfsub$doy),length=100),
                         'year'=rep(2002,100),'bottom_depth'=rep(median(yfsub$bottom_depth,na.rm=TRUE),100),
                         'reg.SST'=mean(yfsub$reg.SST[yfsub$reg.SST<4.175]),'th'="TRUE")
grid.extent3$pred<-predict(thr.pheno,newdata=grid.extent3)
grid.extent3$se<-predict(thr.pheno,newdata=grid.extent3,se=T)[[2]] #select the standard error value from predictions 
grid.extent3$pred.u<-grid.extent3$pred+1.645*grid.extent3$se
grid.extent3$pred.l<-grid.extent3$pred-1.645*grid.extent3$se
grid.extent3$reg.SST<-mean(yfsub$reg.SST[yfsub$reg.SST>4.175])
grid.extent3$th<-"FALSE"
grid.extent3$pred2<-predict(thr.pheno,newdata=grid.extent3)
grid.extent3$se2<-predict(thr.pheno,newdata=grid.extent3,se=T)[[2]]
grid.extent3$pred2.u<-grid.extent3$pred2+1.645*grid.extent3$se2
grid.extent3$pred2.l<-grid.extent3$pred2-1.645*grid.extent3$se2

warmcol<-adjustcolor('orangered3',alpha.f=0.3)
coolcol<-adjustcolor('honeydew2',alpha.f=0.8)

windows(width=12,height=12)
par(mai=c(1,1,0.5,0.5))
plot(grid.extent3$doy,grid.extent3$pred,main='Change in Phenology Due to Two Threshold Conditions',type='l',
     ylim=range(c(grid.extent3$pred.u,grid.extent3$pred2.u,grid.extent3$pred.l,grid.extent3$pred2.l)),
     xlim=range(yfsub$doy),col='black',lwd=2,xlab='Day of the Year',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),cex.lab=1,cex.axis=0.9,cex.main=1)
polygon(c(grid.extent3$doy,rev(grid.extent3$doy)),c(grid.extent3$pred.l,rev(grid.extent3$pred.u)),
        col=coolcol,lty=0)
lines(grid.extent3$doy,grid.extent3$pred,col='grey43',lwd=2)
lines(grid.extent3$doy,grid.extent3$pred2,col='indianred3',lwd=2)
abline(h=0,col='grey79',lty=2,lwd=1.5)
polygon(c(grid.extent3$doy,rev(grid.extent3$doy)),c(grid.extent3$pred2.l,rev(grid.extent3$pred2.u)),
        col=warmcol,lty=0)
legend('bottomright',legend=c(expression(paste(mu, '(<4.175'^0,' C)')),
                              expression(paste(mu,'(>4.175'^0,' C)')),'90% CI','90% CI'),
       col=c('grey43','indianred3','honeydew2',warmcol),pch=c(NA,NA,15,15),lwd=c(2,2,NA,NA),cex=0.8) 


# Predicting Larval Biogeography  -----------------------------------------
#attempting to use above code to predict larval biogeography based on the 2D temp,sal model: 
#using 2006 as the focal year because it has moderate + catches and is ~ in the middle of the time series

#base biogeography: 
nlat=120
nlon=120
latd=seq(min(yflarv.ctd$lat,na.rm=TRUE),max(yflarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(yflarv.ctd$lon,na.rm=TRUE),max(yflarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          yflarv.ctd$lat,yflarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2008)
grid.extent$doy<-as.numeric(median(yflarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(yflarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(yflarv.ctd$lon,na.rm=TRUE),ylim=range(yflarv.ctd$lat,na.rm=TRUE),main='Yellowfin Sole Distribution, Larvae',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(yflarv.ctd$lon[yflarv.ctd$Cper10m2==0],yflarv.ctd$lat[yflarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(yflarv.ctd$lon[yflarv.ctd$Cper10m2>0],
        yflarv.ctd$lat[yflarv.ctd$Cper10m2>0],
        circles=log(yflarv.ctd$Cper10m2+1)[yflarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Improved distribution with temperature,salinity 2D model: 
#plan is to calculate the significant differences in larval biogeography from base model to 2D model
#to show how inclusion of temperature and salinity improves biogeography understanding
nlat=120
nlon=120
latd=seq(min(yflarv.ctd$lat,na.rm=TRUE),max(yflarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(yflarv.ctd$lon,na.rm=TRUE),max(yflarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          yflarv.ctd$lat,yflarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2008)
grid.extent$doy<-as.numeric(median(yflarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(yflarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$se<-predict(lv.base,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA
grid.extent$temperature<-as.numeric(mean(yflarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(yflarv.ctd$salinity))
grid.extent$pred2<-predict(lv.2d,newdata=grid.extent) 
grid.extent$se2<-predict(lv.2d,newdata=grid.extent,se=T)[[2]]
grid.extent$pred2.u<-grid.extent$pred2+1.96*grid.extent$se2
grid.extent$pred2.l<-grid.extent$pred2-1.96*grid.extent$se2
grid.extent$diff<-grid.extent$pred2-grid.extent$pred #calculate differences between base and 2D

grid.extent$sig.pos<-c(grid.extent$pred2.l>grid.extent$pred.u) #isolate areas where there is a higher predicted CPUE in 2D model 
grid.extent$sig.neg<-c(grid.extent$pred2.u<grid.extent$pred.l)
grid.extent$pos.diff<-grid.extent$diff*grid.extent$sig.pos #calculate areas with a significant positive difference in 2D
grid.extent$neg.diff<-grid.extent$diff*grid.extent$sig.neg
max.slope<-max(grid.extent$diff,na.rm=T)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$diff,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(yflarv.ctd$lon,na.rm=TRUE),ylim=range(yflarv.ctd$lat,na.rm=TRUE),
           main=expression(paste('Yellowfin Sole ',Delta,'Larval Distribution w Temp and Salinity')),
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2.8,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#More Simplistic Predicted Larval Biogeography - just plot based on what the 2D model predicts: 
nlat=120
nlon=120
latd=seq(min(yflarv.ctd$lat,na.rm=TRUE),max(yflarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(yflarv.ctd$lon,na.rm=TRUE),max(yflarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          yflarv.ctd$lat,yflarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2008)
grid.extent$doy<-as.numeric(median(yflarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(yflarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(yflarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(yflarv.ctd$salinity))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(yflarv.ctd$lon,na.rm=TRUE),ylim=range(yflarv.ctd$lat,na.rm=TRUE),
           main='Predicted Larval Biogeography, 2D Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Larval Catch Predictions on a Temperature-Salinity Diagram: 
#basically applying same strategy, but instead of a long-lat grid, making a temp-sal grid
ntemp<-100
nsal<-100
tempd<-seq(min(yflarv.ctd$temperature,na.rm=TRUE),max(yflarv.ctd$temperature,na.rm=TRUE),length.out=ntemp)
sald<-seq(min(yflarv.ctd$salinity,na.rm=T),max(yflarv.ctd$salinity,na.rm=T),length.out=nsal)

grid.extent<-expand.grid(sald,tempd)
names(grid.extent)<-c('salinity','temperature')

grid.extent$dist.sal<-NA
grid.extent$dist.temp<-NA
for(k in 1:nrow(grid.extent)){
  dist.sal<-euclidean.distance(grid.extent$salinity[k],yflarv.ctd$salinity[k])
  dist.temp<-euclidean.distance(grid.extent$temperature[k],yflarv.ctd$temperature[k])
  
  grid.extent$dist.sal[k]<-min(dist.sal)
  grid.extent$dist.temp[k]<-min(dist.temp)
}

grid.extent$year<-as.numeric(2005)
grid.extent$lon<-as.numeric(median(yflarv.ctd$lon))
grid.extent$lat<-as.numeric(median(yflarv.ctd$lat))
grid.extent$doy<-as.numeric(median(yflarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(yflarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist.sal>75]<-NA
grid.extent$pred[grid.extent$dist.temp>75]<-NA #check this

windows(width=15,height=15)
par(mai=c(1,1,0.5,0.9))
image.plot(sald,tempd,t(matrix(grid.extent$pred,nrow=length(tempd),ncol=length(sald),byrow=T)),
           col=hcl.colors(100,"PRGn"),xlab='Salinity (psu)',
           ylab=expression(paste("Temperature ("^0, 'C)')),
           xlim=range(yflarv.ctd$salinity,na.rm=T),ylim=range(yflarv.ctd$temperature,na.rm=T),
           main='Larval Biogeography By Temperature and Salinity',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)

