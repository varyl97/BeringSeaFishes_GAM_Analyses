##GAM Predictions 
#visualize results of best model by predicting on a grid: 
#this attempt did not work - 10/11/2021
windows(width=12,height=3.5)
par(mfrow=c(1,4),mai=c(0.7,0.6,0.4,0.4))

#prediction grid
nlat=80
nlon=120
latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calc distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhsub$lat,fhsub$lon)
  grid.extent$dist[k]<-min(dist)
}

#test:     #grid.extent frame needs same variables as model
grid.extent$year<-as.factor(2013) #(below thresh of 1.589 degC)
grid.extent$doy<-as.character(median(fhsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<1.589],na.rm=TRUE)  
grid.extent$temps.in=1.589 #define temps.in = to threshold estimation 
grid.extent$pred<-predict(thr.geo,newdata=grid.extent,
                          se.fit=TRUE,type='response')#threshold geography model 
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhsub$lat,fhsub$lon)
  grid.extent$dist[k]<-min(dist)
}
grid.extent$pred[grid.extent$dist>30000]<-NA

windows()
par(mai=c(0.7,0.6,0.4,0.4))
image(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                         ncol=length(lond),byrow=T)),col=tim.colors(100),ylab="",
      xlab="",xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Distribution',
      cex.main=1.5,cex.lab=1.4,cex.axis=1.4)
symbols(fhsub$lon[fhsub$Cper10m2>0],
        fhsub$lat[fhsub$Cper10m2>0],
        circles=log(fhsub$Cper10m2+1)[fhsub$Cper10m2>0],
        inches=0.1,bg=alpha('grey',f=0.02),fg=alpha('black',f=0.02),add=T)
map("worldHires",fill=T,col="wheat4",add=T)

#Working with Code from Rebecca - 10/11/2021
distances<-function(fhsub,year){
  nlat=80
  nlon=120
  latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
  lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)
  sub_below<-expand.grid(lond,latd)
  names(sub_below)<-c("lon","lat")
  sub_below$year<-2013 #test year
  sub_below$th<-"below"
  sub_below$dist<-NA
  sub_below$doy<-as.character(median(fhsub$doy,na.rm=TRUE))
  sub_below$bottom_depth<-as.character(median(fhsub$bottom_depth,na.rm=TRUE))
  sub_below$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<1.589],na.rm=TRUE)
  sub_below$temps.in<-1.589
  for(i in 1:nrow(sub_below)){
    dist<-distance_function(
      sub_below$lat[i],sub_below$lon[i],
      fhsub$lat,fhsub$lon)
    sub_below$dist[i]<-min(dist)
  }
  return(sub_below)
}

thr.geo<-thr.geo[[best.index.geo]] #function doesn't like the [[best.index]]
tgam_pred<-function(thr.geo,sub_below,year){
  below_pred<-predict(thr.geo,newdata=sub_below,
                      se.fit=TRUE,
                      type='response')
  nlat=80
  nlon=120
  latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
  lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)
  pred_mean_below<-below_pred[[1]]
  pred_se_below<-below_pred[[2]]
  pred_mean_below[sub_below$dist>30000]<-NA
  pred_se_below[sub_below$dist>30000]<-NA
  sub_above<-expand.grid(lond,latd)
  names(sub_above)<-c("lon","lat")
  sub_above$year<-2013
  sub_above$dist<-NA
  sub_above$th<-"above"
  sub_above$doy<-as.character(median(fhsub$doy,na.rm=TRUE))
  sub_above$bottom_depth<-as.character(median(fhsub$bottom_depth,na.rm=TRUE))
  sub_above$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<1.589],na.rm=TRUE)
  sub_above$temps.in<-1.589
  above_pred<-predict(thr.geo,newdata=sub_above,
                      se.fit=TRUE,type='response')
  pred_mean_above<-above_pred[[1]]
  pred_se_above<-above_pred[[2]]
  pred_mean_above[sub_below$dist > 30000] <- NA
  pred_se_above[subset_below$dist > 30000] <- NA
}

#variable-coefficient spatial: 
nlat=80
nlon=120
latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)

grid_vc<-expand.grid(lond,latd)
names(grid_vc)<-c('lon','lat')
grid_vc$dist<-NA
for(k in 1:nrow(grid_vc)){
  dist<-distance.function(grid_vc$lon[k],
                          grid_vc$lat[k],
                          fhsub$lon,
                          fhsub$lat)
  grid_vc$dist[k]<-min(dist)
}

grid_vc$year<-2013
grid_vc$doy<-median(fhsub$doy,na.rm=T)
grid_vc$reg.SST<-mean(fhsub$reg.SST,na.rm=TRUE)
grid_vc$bottom_depth<-as.factor(median(fhsub$bottom_depth,na.rm=T))
grid_vc$pred<-predict(vc.geo,newdata=grid_vc,type='response')
grid_vc$se<-predict(vc.geo,newdata=grid_vc,se=2)[[2]]
grid_vc$pred_up<-grid_vc$pred+1.96*grid_vc$se
grid_vc$pred_lw<-grid_vc$pred-1.96*grid_vc$se
grid_vc$pred[grid_vc$dist>30000]<-NA


windows()
par(mai=c(0.7,0.6,0.4,0.4))
image(lond,latd,t(matrix(grid_vc$pred,nrow=length(latd),
                         ncol=length(lond),byrow=T)),col=tim.colors(100),ylab="",
      xlab="",xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Distribution',
      cex.main=1.5,cex.lab=1.4,cex.axis=1.4)
symbols(fhsub$lon[fhsub$Cper10m2>0],
        fhsub$lat[fhsub$Cper10m2>0],
        circles=log(fhsub$Cper10m2+1)[fhsub$Cper10m2>0],
        inches=0.1,bg=alpha('grey',f=0.02),fg=alpha('black',f=0.02),add=T)
map("worldHires",fill=T,col="wheat4",add=T)























