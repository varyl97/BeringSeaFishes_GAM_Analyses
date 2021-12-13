##GAM Predictions: Northern Rock Sole
#the following code produces predictions and predictive figures based off the best model 
#for Northern Rock Sole, only larvae were evaluated because NRS spawn near/under ice and so their 
#       eggs are rarely collected by ecoFOCI trawls. 

# Preliminary Data and GAM Loading ---------------------------------
nrslarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv',header=TRUE,check.names=TRUE)


lv.base<-readRDS("./GAM Models/nrs_larvae_base.rds")
lv.add.sal<-readRDS("./GAM Models/nrs_larvae_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/nrs_larvae_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/nrs_larvae_addtempsal.rds")
lv.2d<-readRDS("./GAM Models/nrs_larvae_2d.rds")

#get map 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) 

# Predicting Larval Biogeography  -----------------------------------------
#attempting to use above code to predict larval biogeography based on the 2D temp,sal model: 
#using 2005 as the focal year because it has moderate + catches and is ~ in the middle of the time series

#base biogeography: 
nlat=120
nlon=120
latd=seq(min(nrslarv.ctd$lat,na.rm=TRUE),max(nrslarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(nrslarv.ctd$lon,na.rm=TRUE),max(nrslarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          nrslarv.ctd$lat,nrslarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2005)
grid.extent$doy<-as.numeric(median(nrslarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(nrslarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(nrslarv.ctd$lon,na.rm=TRUE),ylim=range(nrslarv.ctd$lat,na.rm=TRUE),main='Northern Rock Sole Distribution, Larvae',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2.4,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(nrslarv.ctd$lon[nrslarv.ctd$Cper10m2==0],nrslarv.ctd$lat[nrslarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(nrslarv.ctd$lon[nrslarv.ctd$Cper10m2>0],
        nrslarv.ctd$lat[nrslarv.ctd$Cper10m2>0],
        circles=log(nrslarv.ctd$Cper10m2+1)[nrslarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Improved distribution with temperature,salinity 2D model: 
#plan is to calculate the significant differences in larval biogeography from base model to 2D model
#to show how inclusion of temperature and salinity improves biogeography understanding
nlat=120
nlon=120
latd=seq(min(nrslarv.ctd$lat,na.rm=TRUE),max(nrslarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(nrslarv.ctd$lon,na.rm=TRUE),max(nrslarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          nrslarv.ctd$lat,nrslarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2005)
grid.extent$doy<-as.numeric(median(nrslarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(nrslarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$se<-predict(lv.base,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA
grid.extent$temperature<-as.numeric(mean(nrslarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(nrslarv.ctd$salinity))
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
           xlim=range(nrslarv.ctd$lon,na.rm=TRUE),ylim=range(nrslarv.ctd$lat,na.rm=TRUE),
           main=expression(paste('Northern Rock Sole ',Delta,'Larval Distribution w Temp and Salinity')),
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2.8,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#More Simplistic Predicted Larval Biogeography - just plot based on what the 2D model predicts: 
nlat=120
nlon=120
latd=seq(min(nrslarv.ctd$lat,na.rm=TRUE),max(nrslarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(nrslarv.ctd$lon,na.rm=TRUE),max(nrslarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          nrslarv.ctd$lat,nrslarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2005)
grid.extent$doy<-as.numeric(median(nrslarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(nrslarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(nrslarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(nrslarv.ctd$salinity))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(nrslarv.ctd$lon,na.rm=TRUE),ylim=range(nrslarv.ctd$lat,na.rm=TRUE),
           main='Predicted Larval Biogeography, 2D Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2.5,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Larval Catch Predictions on a Temperature-Salinity Diagram: 
#basically applying same strategy, but instead of a long-lat grid, making a temp-sal grid
ntemp<-100
nsal<-100
tempd<-seq(min(nrslarv.ctd$temperature,na.rm=TRUE),max(nrslarv.ctd$temperature,na.rm=TRUE),length.out=ntemp)
sald<-seq(min(nrslarv.ctd$salinity,na.rm=T),max(nrslarv.ctd$salinity,na.rm=T),length.out=nsal)

grid.extent<-expand.grid(sald,tempd)
names(grid.extent)<-c('salinity','temperature')

grid.extent$dist.sal<-NA
grid.extent$dist.temp<-NA
for(k in 1:nrow(grid.extent)){
  dist.sal<-euclidean.distance(grid.extent$salinity[k],nrslarv.ctd$salinity[k])
  dist.temp<-euclidean.distance(grid.extent$temperature[k],nrslarv.ctd$temperature[k])
  
  grid.extent$dist.sal[k]<-min(dist.sal)
  grid.extent$dist.temp[k]<-min(dist.temp)
}

grid.extent$year<-as.numeric(2005)
grid.extent$lon<-as.numeric(median(nrslarv.ctd$lon))
grid.extent$lat<-as.numeric(median(nrslarv.ctd$lat))
grid.extent$doy<-as.numeric(median(nrslarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(nrslarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist.sal>75]<-NA
grid.extent$pred[grid.extent$dist.temp>75]<-NA #check this

windows(width=15,height=15)
par(mai=c(1,1,0.5,0.9))
image.plot(sald,tempd,t(matrix(grid.extent$pred,nrow=length(tempd),ncol=length(sald),byrow=T)),
           col=hcl.colors(100,"PRGn"),xlab='Salinity (psu)',
           ylab=expression(paste("Temperature ("^0, 'C)')),
           xlim=range(nrslarv.ctd$salinity,na.rm=T),ylim=range(nrslarv.ctd$temperature,na.rm=T),
           main='Larval Biogeography By Temperature and Salinity',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)

