##GAM Predictions: Flathead Sole
#the following code produces predictions and predictive figures based off the best model 
#for flathead sole, the best model is the threshold geography model wherein spawning geography 
    #varies differently above and below 2.285*C. The best phenology model was also the threshold formulation, 
    #where the threshold temperature was 2.108 *C.  

# Preliminary Data Loading and Formatting ---------------------------------

#link fhsub to regional indices:
fhsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhEggs.csv',header=TRUE,check.names=TRUE) #for egg GAMs
fhlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE) #for larval GAMs

reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(fhsub)){
  fhsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==fhsub$year[i]]}

#load GAMs
eg.base<-readRDS("./GAM Models/fh_egg_base.rds")
thr.pheno<-readRDS("./GAM Models/fh_egg_thr_pheno.rds")
temps.in<-readRDS("./GAM Models/fh_egg_temps_in_pheno.rds") 
best.index.phe<-readRDS("./GAM Models/fh_egg_best_index_pheno.rds")
thr.geo<-readRDS("./GAM Models/fh_egg_thr_geo.rds")
best.index.geo<-readRDS("./GAM Models/fh_egg_best_index_geo.rds")
vc.pheno<-readRDS("./GAM Models/fh_egg_vc_pheno.rds")
vc.geo<-readRDS("./GAM Models/fh_egg_vc_geo.rds") #only using eg.base, thr.pheno, and thr.geo below. 

lv.base<-readRDS("./GAM Models/fh_larv_base.rds")
lv.add.sal<-readRDS("./GAM Models/fh_larv_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/fh_larv_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/fh_larv_tempsal.rds")
lv.2d<-readRDS("./GAM Models/fh_larv_2d.rds") #only using lv.base and lv.2d in this code. 

#get map 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) #seems 'bathy' item is best to add to existing plots
#bathy.dat<-as.bathy(bathy)
#to plot: 
windows()
plot.bathy(bathy,image=T,land=F,n=5,drawlabels=T)


# Visualize on a regular grid ---------------------------------------------
windowsFonts(A="Times New Roman") #for axes labels and plot titles 

#visualize results of best model by predicting on a regular spaced grid: 
#prediction grid
nlat=120
nlon=120
latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhsub$lat,fhsub$lon)
  grid.extent$dist[k]<-min(dist)
}

# Plot Base Geography Model -----------------------------------------------
grid.extent$year<-as.numeric(2013) #for base, just pick a year that has coverage
grid.extent$doy<-as.numeric(median(fhsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(fhsub$reg.SST,na.rm=TRUE) 
grid.extent$pred<-predict(eg.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                         ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
      xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Flathead Sole Distribution, Eggs',
      cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
      legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(fhsub$lon[fhsub$Cper10m2==0],fhsub$lat[fhsub$Cper10m2==0],pch='+',col='white')
symbols(fhsub$lon[fhsub$Cper10m2>0],
        fhsub$lat[fhsub$Cper10m2>0],
        circles=log(fhsub$Cper10m2+1)[fhsub$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)


# Plot Base Phenology Model -----------------------------------------------
grid.extent2<-data.frame('lon'=rep(-155,100),
                         'lat'=rep(51,100),'doy'=seq(min(fhsub$doy),max(fhsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(fhsub$bottom_depth,na.rm=TRUE),100))#setting up another clean grid on which to predict
grid.extent2$pred<-predict(eg.base,newdata=grid.extent2)
grid.extent2$se<-predict(eg.base,newdata=grid.extent2,se=T)[[2]] #select the standard error 
grid.extent2$pred.u<-grid.extent2$pred+(1.96*grid.extent2$se) #calculate upper confidence interval (95% CI)
grid.extent2$pred.l<-grid.extent2$pred-(1.96*grid.extent2$se) #calculate lower CI 


windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
plot(grid.extent2$doy,grid.extent2$pred,main='Flathead Sole Base Phenology, Eggs',type='l',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),xlab='Day of Year',cex.lab=1,
     cex.axis=1,cex.main=1,cex.axis=0.9,xlim=range(fhsub$doy),
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred.l)),
     col='grey1',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col='azure3',lty=0)
lines(grid.extent2$doy,grid.extent2$pred,col='grey1')
legend('topleft',legend=c(expression(paste("(log(C/(10m"^2,')+1)')),'95% CI'),
       col=c('grey1','azure3'),pch=c(NA,15),lty=c(1,NA),lwd=2,cex=1)
abline(h=0,col='grey79',lty=2,lwd=1.5)

#using base graphics: 
windows()
par(mai=c(1,1,0.5,0.5))
plot(eg.base,select=2,main='Flathead Sole Base Phenology, Eggs',
     seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies (edf: 8.79)',ylim=c(-3,3))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)


#TEMP EFFECT: Calculate Differences Due to Different Temperature Regimes Based on Best Model --------
#start with threshold geography model to find differences between two predictions to calculate local slopes 
nlat=120
nlon=120
latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhsub$lat,fhsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-2013
grid.extent$doy<-median(fhsub$doy)
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<temps.in[best.index.geo]]) #threshold temp chosen by AIC values
grid.extent$th<-"TRUE"
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=T)
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)
grid.extent$se<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se #95% CI here
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA #remove predictions that are too far from positive data values
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>temps.in[best.index.geo]])
grid.extent$th<-"FALSE"
grid.extent$pred2<-predict(thr.geo,newdata=grid.extent)
grid.extent$se2<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
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
           xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Change in FHS(E) Distribution w Threshold Temperature Effect',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),
           legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)#would prefer to have legend within plot margins, and for all font to be times, but not sure how to do that. 
map("worldHires",fill=T,col="seashell2",add=T)

#now add in the Phenology effect from this model (again, using this model because it produced most deviance explained and lowest AIC): 
#using base graphics here, no need to overlay anything
windows()
par(mai=c(1,1,0.5,0.5))
plot(thr.geo,select=1,main='Flathead Sole Threshold Geo Phenology, Eggs',
     seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies (edf: 8.78)',ylim=c(-1.5,3))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)


#plot the two phenology smooths together, one from the base model and one from the threshold geography model to see the temp effect: 
col<-adjustcolor('tomato4',alpha.f=0.3)

windows()
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5))
plot(eg.base,select=2,main='Flathead Sole Phenology, Eggs',seWithMean=TRUE,
     ylim=c(-1.5,3))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5),new=TRUE)
plot(thr.geo,select=1,seWithMean=TRUE,shade=TRUE,shade.col=col,ylim=c(-1.5,3))
legend('topright',legend=c('Base','Threshold Geography'),col=c(NA,col),lwd=c(2,2),cex=0.8)
mtext(c("Day of Year","Anomalies in log(CPUE+1)"),side=c(1,2),line=2.5)

#Plot Threshold Phenology (best model after new iteration, 12/13/2021):
col<-adjustcolor('tomato4',alpha.f=0.3)

windows()
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5))
plot(thr.pheno,select=4,main='Flathead Sole Phenology, Threshold Effect',seWithMean=TRUE,
     ylim=c(-3,3))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5),new=TRUE)
plot(thr.pheno,select=3,seWithMean=TRUE,shade=TRUE,shade.col=col,ylim=c(-3,3))
legend('topright',legend=c('Below 2.12*C','Above 2.12*C'),col=c('grey25',col),lwd=c(2,2),cex=0.8)
mtext(c("Day of Year","Anomalies in log(CPUE+1)"),side=c(1,2),line=2.5)


# Plot Geography from Phenology Model -----------------------------------------------
grid.extent$year<-as.numeric(2013) #for base, just pick a year that has coverage
grid.extent$doy<-as.numeric(median(fhsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(fhsub$reg.SST,na.rm=TRUE) 
grid.extent$th<-as.character("TRUE")
grid.extent$pred<-predict(thr.pheno,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Flathead Sole Distribution, Th.Ph, Eggs',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#RELIC: For added information: threshold phenology prediction (this was the best model to explain phenology; In new iteration (12/7/2021, this is best model)
grid.extent3<-data.frame('lon'=rep(-155,100),'lat'=rep(51,100),'doy'=seq(min(fhsub$doy),max(fhsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(fhsub$bottom_depth,na.rm=TRUE),100),
                         'reg.SST'=mean(fhsub$reg.SST[fhsub$reg.SST<2.108]),'th'="TRUE")
grid.extent3$pred<-predict(thr.pheno,newdata=grid.extent3)
grid.extent3$se<-predict(thr.pheno,newdata=grid.extent3,se=T)[[2]] #select the standard error value from predictions 
grid.extent3$pred.u<-grid.extent3$pred+1.645*grid.extent3$se
grid.extent3$pred.l<-grid.extent3$pred-1.645*grid.extent3$se
grid.extent3$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>2.108])
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
     xlim=range(fhsub$doy),col='black',lwd=2,xlab='Day of the Year',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),cex.lab=1,cex.axis=0.9,cex.main=1)
polygon(c(grid.extent3$doy,rev(grid.extent3$doy)),c(grid.extent3$pred.l,rev(grid.extent3$pred.u)),
        col=coolcol,lty=0)
lines(grid.extent3$doy,grid.extent3$pred,col='grey43',lwd=2)
lines(grid.extent3$doy,grid.extent3$pred2,col='indianred3',lwd=2)
abline(h=0,col='grey79',lty=2,lwd=1.5)
polygon(c(grid.extent3$doy,rev(grid.extent3$doy)),c(grid.extent3$pred2.l,rev(grid.extent3$pred2.u)),
        col=warmcol,lty=0)
legend('topleft',legend=c(expression(paste(mu, '(<2.108'^0,' C)')),
                          expression(paste(mu,'(>2.108'^0,' C)')),'90% CI','90% CI'),
       col=c('grey43','indianred3','honeydew2',warmcol),pch=c(NA,NA,15,15),lwd=c(2,2,NA,NA),cex=0.8) 


# Predicting Larval Biogeography  -----------------------------------------
#attempting to use above code to predict larval biogeography based on the 2D temp,sal model: 
#using 2006 as the focal year because it has moderate + catches and is ~ in the middle of the time series

#base biogeography: 
nlat=120
nlon=120
latd=seq(min(fhlarv.ctd$lat,na.rm=TRUE),max(fhlarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(fhlarv.ctd$lon,na.rm=TRUE),max(fhlarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhlarv.ctd$lat,fhlarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(fhlarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(fhlarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(fhlarv.ctd$lon,na.rm=TRUE),ylim=range(fhlarv.ctd$lat,na.rm=TRUE),main='Flathead Sole Distribution, Larvae',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(fhlarv.ctd$lon[fhlarv.ctd$Cper10m2==0],fhlarv.ctd$lat[fhlarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(fhlarv.ctd$lon[fhlarv.ctd$Cper10m2>0],
        fhlarv.ctd$lat[fhlarv.ctd$Cper10m2>0],
        circles=log(fhlarv.ctd$Cper10m2+1)[fhlarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Improved Larval Biogeography: based on what the 2D model predicts: 
nlat=120
nlon=120
latd=seq(min(fhlarv.ctd$lat,na.rm=TRUE),max(fhlarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(fhlarv.ctd$lon,na.rm=TRUE),max(fhlarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhlarv.ctd$lat,fhlarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(fhlarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(fhlarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(fhlarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(fhlarv.ctd$salinity))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(fhlarv.ctd$lon,na.rm=TRUE),ylim=range(fhlarv.ctd$lat,na.rm=TRUE),
           main='Predicted Larval Biogeography, 2D Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Larval Catch Predictions on a Temperature-Salinity Diagram: 
#basically applying same strategy, but instead of a long-lat grid, making a temp-sal grid
ntemp<-100
nsal<-100
tempd<-seq(min(fhlarv.ctd$temperature,na.rm=TRUE),max(fhlarv.ctd$temperature,na.rm=TRUE),length.out=ntemp)
sald<-seq(min(fhlarv.ctd$salinity,na.rm=T),max(fhlarv.ctd$salinity,na.rm=T),length.out=nsal)

grid.extent<-expand.grid(sald,tempd)
names(grid.extent)<-c('salinity','temperature')

grid.extent$dist.sal<-NA
grid.extent$dist.temp<-NA
for(k in 1:nrow(grid.extent)){
  dist.sal<-euclidean.distance(grid.extent$salinity[k],
                               fhlarv.ctd$salinity[fhlarv.ctd$Cper10m2>0][k])
  dist.temp<-euclidean.distance(grid.extent$temperature[k],
                                fhlarv.ctd$temperature[fhlarv.ctd$Cper10m2>0][k])
  
  grid.extent$dist.sal[k]<-min(dist.sal)
  grid.extent$dist.temp[k]<-min(dist.temp)
}

grid.extent$year<-as.numeric(2005)
grid.extent$lon<-as.numeric(median(fhlarv.ctd$lon))
grid.extent$lat<-as.numeric(median(fhlarv.ctd$lat))
grid.extent$doy<-as.numeric(median(fhlarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(fhlarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist.sal>1.13]<-NA
grid.extent$pred[grid.extent$dist.temp>6.86]<-NA #threshold based on 3rd quartile

windows(width=15,height=15)
par(mai=c(1,1,0.5,0.9))
image.plot(sald,tempd,t(matrix(grid.extent$pred,nrow=length(tempd),ncol=length(sald),byrow=T)),
           col=hcl.colors(100,"PRGn"),xlab='Salinity (psu)',
           ylab=expression(paste("Temperature ("^0, 'C)')),
           xlim=range(fhlarv.ctd$salinity,na.rm=T),ylim=range(fhlarv.ctd$temperature,na.rm=T),
           main='Larval Biogeography By Temperature and Salinity',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
symbols(fhlarv.ctd$salinity[fhlarv.ctd$Cper10m2>0],
        fhlarv.ctd$temperature[fhlarv.ctd$Cper10m2>0],
        circles=log(fhlarv.ctd$Cper10m2+1)[fhlarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)

















