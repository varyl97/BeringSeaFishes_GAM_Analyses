# Generate Presence-Absence Models for Larval Biogeography ----------------
#the following code will generate presence (1) or absence (0) indices for all lines of larval data
#then, a generalized additive model is fit using a binomial distribution 
#this code includes both straight model outputs and predictions based on model (to improve figure visualization)

str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) 

# Applying Presence-Absence Values to Data --------------------------------

aplarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE)
fhlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE)
nrslarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pclarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pklarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',header=TRUE,check.names=TRUE)
rxlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_RxLarv_wCTD.csv',header=TRUE,check.names=TRUE)
yflarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv',header=TRUE,check.names=TRUE)

aplarv.ctd <- aplarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                    Cper10m2<=0 ~ as.numeric('0')))
fhlarv.ctd <- fhlarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2<=0 ~ as.numeric('0')))
nrslarv.ctd <- nrslarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                    Cper10m2<=0 ~ as.numeric('0')))
pclarv.ctd <- pclarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2<=0 ~ as.numeric('0')))
pklarv.ctd <- pklarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2<=0 ~ as.numeric('0')))
rxlarv.ctd <- rxlarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2<=0 ~ as.numeric('0')))
yflarv.ctd <- yflarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2<=0 ~ as.numeric('0')))#making a new variable, "obs", that is binary for presence (1) or absence (0)


# Generalized Additive Model Formulations ---------------------------------

## Flathead Sole, 
#models: 
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=fhlarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=fhlarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=fhlarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=fhlarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.fhlarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.fhlarv$AIC_value,main='AIC Results for fh Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.fhlarv$AIC_value,labels=round(aic.fhlarv$AIC_value),pos=c(4,3,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with the best model: (2d for FHS)
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
grid.extent$salinity<-as.numeric(mean(fhlarv.ctd$salinity,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(fhlarv.ctd$temperature,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(fhlarv.ctd$lon,na.rm=TRUE),ylim=range(fhlarv.ctd$lat,na.rm=TRUE),
           main='Flathead Sole Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

## Rex Sole, which had the worst performance of larval models:  
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=rxlarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=rxlarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=rxlarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=rxlarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.rxlarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.rxlarv$AIC_value,main='AIC Results for rx Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.rxlarv$AIC_value,labels=round(aic.rxlarv$AIC_value),pos=c(4,3,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, 2d: 
nlat=120
nlon=120
latd=seq(min(rxlarv.ctd$lat,na.rm=TRUE),max(rxlarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(rxlarv.ctd$lon,na.rm=TRUE),max(rxlarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          rxlarv.ctd$lat,rxlarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(rxlarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(rxlarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$salinity<-as.numeric(mean(rxlarv.ctd$salinity,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(rxlarv.ctd$temperature,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(rxlarv.ctd$lon,na.rm=TRUE),ylim=range(rxlarv.ctd$lat,na.rm=TRUE),
           main='Rex Sole Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(rxlarv.ctd$lon[rxlarv.ctd$Cper10m2==0],rxlarv.ctd$lat[rxlarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(rxlarv.ctd$lon[rxlarv.ctd$Cper10m2>0],
        rxlarv.ctd$lat[rxlarv.ctd$Cper10m2>0],
        circles=log(rxlarv.ctd$Cper10m2+1)[rxlarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

## Alaska Plaice: 
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=aplarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=aplarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=aplarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=aplarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.aplarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.aplarv$AIC_value,main='AIC Results for ap Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.aplarv$AIC_value,labels=round(aic.aplarv$AIC_value),pos=c(4,4,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, add temp: 
nlat=120
nlon=120
latd=seq(min(aplarv.ctd$lat,na.rm=TRUE),max(aplarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(aplarv.ctd$lon,na.rm=TRUE),max(aplarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          aplarv.ctd$lat,aplarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(aplarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(aplarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(aplarv.ctd$temperature,na.rm=TRUE))
grid.extent$pred<-predict(lv.add.temp,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(aplarv.ctd$lon,na.rm=TRUE),ylim=range(aplarv.ctd$lat,na.rm=TRUE),
           main='Alaska Plaice Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

## Walleye Pollock
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=pklarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=pklarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=pklarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=pklarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.pklarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.pklarv$AIC_value,main='AIC Results for pk Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.pklarv$AIC_value,labels=round(aic.pklarv$AIC_value),pos=c(4,4,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, 2d: 
nlat=120
nlon=120
latd=seq(min(pklarv.ctd$lat,na.rm=TRUE),max(pklarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(pklarv.ctd$lon,na.rm=TRUE),max(pklarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          pklarv.ctd$lat,pklarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(pklarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(pklarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(pklarv.ctd$temperature,na.rm=TRUE))
grid.extent$salinity<-as.numeric(mean(pklarv.ctd$salinity,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(pklarv.ctd$lon,na.rm=TRUE),ylim=range(pklarv.ctd$lat,na.rm=TRUE),
           main='Walleye Pollock Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)

## Yellowfin Sole: 
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=yflarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=yflarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=yflarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=yflarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.yflarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.yflarv$AIC_value,main='AIC Results for yf Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.yflarv$AIC_value,labels=round(aic.yflarv$AIC_value),pos=c(4,4,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, 2d (tied with add temp):
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

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(yflarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(yflarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(yflarv.ctd$temperature,na.rm=TRUE))
grid.extent$salinity<-as.numeric(mean(yflarv.ctd$salinity,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(yflarv.ctd$lon,na.rm=TRUE),ylim=range(yflarv.ctd$lat,na.rm=TRUE),
           main='Yellowfin Sole Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)


## Northern Rock Sole: 
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=nrslarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=nrslarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=nrslarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=nrslarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.nrslarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.2d))

windows()
plot(c(1:4),aic.nrslarv$AIC_value,main='AIC Results for nrs Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.nrslarv$AIC_value,labels=round(aic.nrslarv$AIC_value),pos=c(4,4,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, 2d:
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

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(nrslarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(nrslarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(nrslarv.ctd$temperature,na.rm=TRUE))
grid.extent$salinity<-as.numeric(mean(nrslarv.ctd$salinity,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(nrslarv.ctd$lon,na.rm=TRUE),ylim=range(nrslarv.ctd$lat,na.rm=TRUE),
           main='Northern Rock Sole Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)


## Pacific Cod:
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=pclarv.ctd,family=binomial)
summary(lv.base)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=pclarv.ctd,family=binomial)
summary(lv.add.sal)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=pclarv.ctd,family=binomial)
summary(lv.add.temp)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=pclarv.ctd,family=binomial)
summary(lv.2d)

#AIC checking to identify best model: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.2d<-AIC(lv.2d)

aic.pclarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                  'Sal-Temp 2D'),
                        'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                      aic.2d))

windows()
plot(c(1:4),aic.pclarv$AIC_value,main='AIC Results for pc Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.pclarv$AIC_value,labels=round(aic.pclarv$AIC_value),pos=c(4,4,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF","#FDE725FF"),
       lwd=3,lty=1)

#prediction with best model, 2d:
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

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(pclarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(pclarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$temperature<-as.numeric(mean(pclarv.ctd$temperature,na.rm=TRUE))
grid.extent$salinity<-as.numeric(mean(pclarv.ctd$salinity,na.rm=TRUE))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent,type="response")
grid.extent$pred[grid.extent$dist>30000]<-NA

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(pclarv.ctd$lon,na.rm=TRUE),ylim=range(pclarv.ctd$lat,na.rm=TRUE),
           main='Pacific Cod Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="seashell2",add=T)












