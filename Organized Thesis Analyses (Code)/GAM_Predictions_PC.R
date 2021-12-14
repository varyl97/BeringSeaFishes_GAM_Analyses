##GAM Predictions: Pacific Cod
#the following code produces predictions and predictive figures based off the best model 
#for Pacific cod, only larvae were evaluated because Pacific cod spawn near/under ice and so their 
#       eggs are rarely collected by ecoFOCI trawls. 

# Preliminary Data and GAM Loading ---------------------------------
pclarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv',header=TRUE,check.names=TRUE)

lv.base<-readRDS("./GAM Models/pc_larvae_base.rds")
lv.add.sal<-readRDS("./GAM Models/pc_larvae_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/pc_larvae_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/pc_larvae_addtempsal.rds")
lv.2d<-readRDS("./GAM Models/pc_larvae_2d.rds")

#get map 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) 

# Predicting Larval Biogeography  -----------------------------------------
#attempting to use above code to predict larval biogeography based on the 2D temp,sal model: 
#using 2005 as the focal year because it has moderate + catches and is ~ in the middle of the time series

#base biogeography: 
nlat=120
nlon=120
latd=seq(min(pclarv.ctd$lat,na.rm=TRUE),max(pclarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(pclarv.ctd$lon,na.rm=TRUE),max(pclarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          pclarv.ctd$lat,pclarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2005)
grid.extent$doy<-as.numeric(median(pclarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(pclarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(pclarv.ctd$lon,na.rm=TRUE),ylim=range(pclarv.ctd$lat,na.rm=TRUE),main='Pacific Cod Distribution, Larvae',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(pclarv.ctd$lon[pclarv.ctd$Cper10m2==0],pclarv.ctd$lat[pclarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(pclarv.ctd$lon[pclarv.ctd$Cper10m2>0],
        pclarv.ctd$lat[pclarv.ctd$Cper10m2>0],
        circles=log(pclarv.ctd$Cper10m2+1)[pclarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Predicted Larval Biogeography - based on what the 2D model predicts: 
nlat=120
nlon=120
latd=seq(min(pclarv.ctd$lat,na.rm=TRUE),max(pclarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(pclarv.ctd$lon,na.rm=TRUE),max(pclarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          pclarv.ctd$lat,pclarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2005)
grid.extent$doy<-as.numeric(median(pclarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(pclarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(pclarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(pclarv.ctd$salinity))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(pclarv.ctd$lon,na.rm=TRUE),ylim=range(pclarv.ctd$lat,na.rm=TRUE),
           main='Predicted Larval Biogeography, 2D Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

#Larval Catch Predictions on a Temperature-Salinity Diagram: 
#basically applying same strategy, but instead of a long-lat grid, making a temp-sal grid
ntemp<-100
nsal<-100
tempd<-seq(min(pclarv.ctd$temperature,na.rm=TRUE),max(pclarv.ctd$temperature,na.rm=TRUE),length.out=ntemp)
sald<-seq(min(pclarv.ctd$salinity,na.rm=T),max(pclarv.ctd$salinity,na.rm=T),length.out=nsal)

grid.extent<-expand.grid(sald,tempd)
names(grid.extent)<-c('salinity','temperature')

grid.extent$dist.sal<-NA
grid.extent$dist.temp<-NA
for(k in 1:nrow(grid.extent)){
  dist.sal<-euclidean.distance(grid.extent$salinity[k],
                               pclarv.ctd$salinity[pclarv.ctd$Cper10m2>0][k])
  dist.temp<-euclidean.distance(grid.extent$temperature[k],
                                pclarv.ctd$temperature[pclarv.ctd$Cper10m2>0][k])
  
  grid.extent$dist.sal[k]<-min(dist.sal)
  grid.extent$dist.temp[k]<-min(dist.temp)
}

grid.extent$year<-as.numeric(2005)
grid.extent$lon<-as.numeric(median(pclarv.ctd$lon))
grid.extent$lat<-as.numeric(median(pclarv.ctd$lat))
grid.extent$doy<-as.numeric(median(pclarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(pclarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist.sal>1.917]<-NA
grid.extent$pred[grid.extent$dist.temp>6.975]<-NA #check this

windows(width=15,height=15)
par(mai=c(1,1,0.5,0.9))
image.plot(sald,tempd,t(matrix(grid.extent$pred,nrow=length(tempd),ncol=length(sald),byrow=T)),
           col=hcl.colors(100,"PRGn"),xlab='Salinity (psu)',
           ylab=expression(paste("Temperature ("^0, 'C)')),
           xlim=range(pclarv.ctd$salinity,na.rm=T),ylim=range(pclarv.ctd$temperature,na.rm=T),
           main='Larval Biogeography By Temperature and Salinity',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
symbols(pclarv.ctd$salinity[pclarv.ctd$Cper10m2>0],
        pclarv.ctd$temperature[pclarv.ctd$Cper10m2>0],
        circles=log(pclarv.ctd$Cper10m2+1)[pclarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)

