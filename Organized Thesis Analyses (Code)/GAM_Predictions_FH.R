##GAM Predictions 
#the following code produces predictions and predictive figures based off the best model 
#for flathead sole, the best model is the threshold geography model wherein spawning geography 
    #varies differently above and below 2.285*C. 

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
vc.geo<-readRDS("./GAM Models/fh_egg_vc_geo.rds")

#get map 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name)
bathy.dat<-as.bathy(bathy)
#to plot: 
windows()
plot.bathy(bathy.dat,image=T,land=F,n=5,drawlabels=T) #note: this works as an argument alone, but can't seem to add this to existing plots 

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

windows()
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                         ncol=length(lond),byrow=T)),col=viridis(100),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
      xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Flathead Sole Distribution',
      cex.main=1.5,cex.lab=1.4,cex.axis=1.4,legend.line=2,legend.lab=expression("Predicted Log(Catch Per 10m^2+1)"))
points(fhsub$lon[fhsub$Cper10m2==0],fhsub$lat[fhsub$Cper10m2==0],pch='+',col='white')
symbols(fhsub$lon[fhsub$Cper10m2>0],
        fhsub$lat[fhsub$Cper10m2>0],
        circles=log(fhsub$Cper10m2+1)[fhsub$Cper10m2>0],
        inches=0.1,bg='grey',fg='black',add=T)
map("worldHires",fill=T,col="burlywood1",add=T)


# Plot Base Phenology Model -----------------------------------------------

grid.extent2<-data.frame('lon'=rep(-155,100),
                         'lat'=rep(51,100),'doy'=seq(min(fhsub$doy),max(fhsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(fhsub$bottom_depth,na.rm=TRUE),100))#setting up another clean grid on which to predict
grid.extent2$pred<-predict(eg.base,newdata=grid.extent2)
grid.extent2$se<-predict(eg.base,newdata=grid.extent2,se=T)[[2]] #select the standard error 
grid.extent2$pred.u<-grid.extent2$pred+1.645*grid.extent2$se #calculate upper confidence interval (90% CI)
grid.extent2$pred.l<-grid.extent2$pred-1.645*grid.extent2$se #calculate lower CI 

windows()
par(mai=c(1,1,0.5,0.9))
plot(grid.extent2$doy,grid.extent2$pred,main='Flathead Sole Phenology',type='l',
     ylab=expression('Egg Density (log(Catch per 10m^2+1)'),xlab='Day of Year',cex.lab=1.4,
     cex.axis=1.4,cex.main=1.5,xlim=range(fhsub$doy),
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred.l)),
     col='steelblue2',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col='steelblue2',lty=0)
abline(h=0,col='sienna3',lty=2,lwd=2)

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
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST<2.285]) #threshold temp chosen by AIC values
grid.extent$th<-"TRUE"
grid.extent$bottom_depth<-median(fhsub$bottom_depth,na.rm=T)
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)
grid.extent$se<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se #95% CI here
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA #remove predictions that are too far from positive data values
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>2.285])
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

windows(width=18,height=9)
par(mfrow=c(1,2),mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$diff,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=hcl.colors(100,"PRGn"),ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')), #PRGn diverges more clearly, helping interpretation
           xlim=range(fhsub$lon),ylim=range(fhsub$lat),main='Change in FHS Distribution Due to Temperature',
           cex.main=1,cex.lab=1,cex.axis=0.9,horizontal=1,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),
           legend.shrink=0.3) #would prefer to have legend within plot margins, and for all font to be times, but not sure how to do that. 
map("worldHires",fill=T,col="slategray2",add=T)

#now add in the phenology effect from this model: 
grid.extent2<-data.frame('lon'=rep(-155,100),'lat'=rep(51,100),'doy'=seq(min(fhsub$doy),max(fhsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(fhsub$bottom_depth,na.rm=TRUE),100),
                         'reg.SST'=mean(fhsub$reg.SST[fhsub$reg.SST<2.285]),'th'="TRUE")
grid.extent2$pred<-predict(thr.geo,newdata=grid.extent2)
grid.extent2$se<-predict(thr.geo,newdata=grid.extent2,se=T)[[2]] #select the standard error value from predictions 
grid.extent2$pred.u<-grid.extent2$pred+1.645*grid.extent2$se
grid.extent2$pred.l<-grid.extent2$pred-1.645*grid.extent2$se
grid.extent2$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>2.285])
grid.extent2$th<-"FALSE"
grid.extent2$pred2<-predict(thr.geo,newdata=grid.extent2)
grid.extent2$se2<-predict(thr.geo,newdata=grid.extent2,se=T)[[2]]
grid.extent2$pred2.u<-grid.extent$pred2+1.645*grid.extent2$se2
grid.extent2$pred2.l<-grid.extent$pred2-1.645*grid.extent2$se2

plot(grid.extent2$doy,grid.extent$pred,main='Change in Phenology Due to Two Threshold Conditions',type='l',
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred2.u,grid.extent2$pred.l,grid.extent2$pred2.l)),
     xlim=range(fhsub$doy),col='steelblue2',lwd=2,xlab='Day of the Year',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),cex.lab=1,cex.axis=0.9,cex.main=1)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col='steelblue2',lty=0,alpha=0.5)
lines(grid.extent2$doy,grid.extent2$pred2,col='violetred2',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred2.l,rev(grid.extent2$pred2.u)),
        col='violetred2',lty=0,alpha=0.5)
legend('bottomright',legend=c('\mu (<2.285)','\mu(> 2.285)'))


#now look at changes in the phenology model - using threshold phenology because it has a lower AIC than variable-coefficient phenology model 
grid.extent<-data.frame('lon'=rep(-155,100),'lat'=rep(51,100),'doy'=seq(min(fhsub$doy),max(fhsub$doy),length=100),
                        'reg.SST'=rep(mean(fhsub$reg.SST[fhsub$reg.SST<2.285]),100),
                        'bottom_depth'=rep(median(fhsub$bottom_depth,na.rm=TRUE),100),
                        'th'=rep("TRUE"),100) #create a data.frame with the same variables as those included in the GAM 
grid.extent$pred<-predict(thr.pheno,newdata=grid.extent)
grid.extent$se<-predict(thr.pheno,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.645*grid.extent$se
grid.extent$pred.l<-grid.extent$pred-1.645*grid.extent$se
grid.extent$reg.SST<-mean(fhsub$reg.SST[fhsub$reg.SST>2.285])
grid.extent$th<-"FALSE"
grid.extent$pred2<-predict(thr.pheno,newdata=grid.extent)
grid.extent$se2<-predict(thr.pheno,newdata=grid.extent)[[2]]
grid.extent$pred2.u<-grid.extent$pred2+1.645*grid.extent$se2
grid.extent$pred2.l<-grid.extent$pred2-1.645*grid.extent$se2

plot(grid.extent$doy,grid.extent$pred,main='Change in FHS Phenology Due to Temperature',type='l',
     ylab=expression('Egg Density (log(Catch per 10m^2+1)'),xlab='Day of Year',cex.lab=1,
     cex.axis=0.9,cex.main=1,xlim=range(fhsub$doy),
     ylim=range(c(grid.extent$pred.u,grid.extent$pred2.u,
                  grid.extent$pred.l,grid.extent$pred2.l)),
     xlim=range(fhsub$doy),
     col='steelblue2',lwd=2)
polygon(c(grid.extent$doy,rev(grid.extent$doy)),c(grid.extent$pred.l,rev(grid.extent$pred.u)),
        col='steelblue2',lty=0)
lines(grid.extent$doy,grid.extent$pred2,col='violetred2',lwd=2)
polygon(c(grid.extent$doy,rev(grid.extent$doy)),c(grid.extent$pred2.lw,rev(grid.extent$pred2.u)),
        col='violetred2',lty=0)
legend("bottomright",legend=c('Mean SST Below 2.285','Mean SST Above 2.285'),
       col=c('steelblue2','violetred2'),lty=1,lwd=2,cex=1)
abline(h=0,col='sienna3',lty=2,lwd=2)

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























