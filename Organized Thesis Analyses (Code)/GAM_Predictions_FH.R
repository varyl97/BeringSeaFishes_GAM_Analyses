##GAM Predictions 

# Initial Attempts - Code based off of Lorenzo's --------------------------
#link fhsub to regional indices: 

reg.SST<-read.csv('../Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.SST) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(fhsub)){
  fhsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==fhsub$year[i]]}


#visualize results of best model by predicting on a grid: 
#this attempt did not work - 10/11/2021
windows(width=12,height=3.5)
par(mfrow=c(1,4),mai=c(0.7,0.6,0.4,0.4))

#prediction grid
nlat=120
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
# Threshold Geo Model. Make sure grid.extent has same variables as model --------

grid.extent$year<-as.numeric(2013) #(below thresh of 1.24 degC)
grid.extent$doy<-as.numeric(median(fhsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<1.24],na.rm=TRUE)  
grid.extent$th=as.character("TRUE")
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)





windows()
par(mai=c(0.7,0.6,0.4,0.4))
image(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                         ncol=length(lond),byrow=T)),col=tim.colors(100),ylab="",
      xlab="",xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Distribution',
      cex.main=1.5,cex.lab=1.4,cex.axis=1.4)
symbols(fhsub$lon[fhsub$Cper10m2>0&fhsub$year==2013],#this isn't working for some reason. 
        fhsub$lat[fhsub$Cper10m2>0&fhsub$year==2013],
        circles=log(fhsub$Cper10m2+1)[fhsub$Cper10m2>0&fhsub$year==2013],
        inches=0.1,bg='grey',fg='black',add=T)
map("worldHires",fill=T,col="wheat4",add=T)

#trying above code with the base model: 
eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=fhsub,family=tw(link='log'),method='REML')

nlat=80
nlon=120
latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          fhsub$lat,fhsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.factor(2013) #(below thresh of 1.24 degC)
grid.extent$doy<-as.character(median(fhsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=TRUE)
grid.extent$pred<-predict(eg.base,newdata=grid.extent,
                          se.fit=TRUE,type='response')



# Working with Code from Rebecca - 10/11/2021 ----------------------------
#write a function instead 
distances<-function(fhsub,year){
  nlat=80
  nlon=120
  latd=seq(min(fhsub$lat),max(fhsub$lat),length.out=nlat)
  lond=seq(min(fhsub$lon),max(fhsub$lon),length.out=nlon)
  sub_below<-expand.grid(lond,latd) #"below" = below threshold temp 
  names(sub_below)<-c("lon","lat")
  sub_below$year<-2013 #test year
  sub_below$th<-as.factor("TRUE")
  sub_below$dist<-NA
  sub_below$doy<-as.character(median(fhsub$doy,na.rm=TRUE))
  sub_below$bottom_depth<-as.character(median(fhsub$bottom_depth,na.rm=TRUE))
  sub_below$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<1.589],na.rm=TRUE)
  for(i in 1:nrow(sub_below)){
    dist<-distance_function(
      sub_below$lat[i],sub_below$lon[i],
      fhsub$lat,fhsub$lon)
    sub_below$dist[i]<-min(dist)
  }
  return(sub_below) #when I try head(sub_below), object is not found 
}

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
  sub_above$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>1.589],na.rm=TRUE)
  above_pred<-predict(thr.geo,newdata=sub_above,
                      se.fit=TRUE,type='response')
  pred_mean_above<-above_pred[[1]]
  pred_se_above<-above_pred[[2]]
  pred_mean_above[sub_below$dist > 30000] <- NA
  pred_se_above[subset_below$dist > 30000] <- NA
  pred_mean_up_below<-pred_mean_below+1.96*pred_se_below #calculate 95% CI for each prediction 
  pred_mean_down_below<-pred_mean_below-1.96*pred_se_below
  pred_mean_up_above<-pred_mean_above+1.96*pred_se_below
  pred_mean_down_above<-pred_mean_above-1.96*pred_se_below
  significant_low<-pred_mean_up_above<pred_mean_down_below
  significant_high<-pred_mean_down_above>pred_mean_up_below #significant differences in preds 
  return(list(significant_high,significant_low,before_prediction,after_prediction))
} #also not sure this is working. I think this above code may just be above my paygrade at this point.. I can't easily isolate what 
#   within this is causing errors but the above functions do not result in anything that I can figure out how to plot. 

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
#results in this error: 'Error in Predict.matrix.tprs.smooth(object, dk$data) : 
#NA/NaN/Inf in foreign function call (arg 1)
#In addition: Warning message:
  #In Ops.factor(xx, object$shift[i]) : '-' not meaningful for factors'
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























