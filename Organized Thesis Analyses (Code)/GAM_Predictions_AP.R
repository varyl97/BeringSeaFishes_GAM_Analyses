##GAM Predictions: Alaska Plaice
#the following code produces predictions and predictive figures based off the best model 
#for Alaska plaice, the best model is the threshold geography model wherein spawning geography 
#varies differently above and below 2.064 *C. 

# Preliminary Data Loading and Formatting ---------------------------------

#link apsub to regional indices:
apsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApEggs.csv',header=TRUE,check.names=TRUE) #for egg GAMs
aplarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE) #for larval GAMs

reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(apsub)){
  apsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==apsub$year[i]]}

#load GAMs
eg.base<-readRDS("./GAM Models/ap_egg_base.rds")
thr.pheno<-readRDS("./GAM Models/ap_egg_thr_pheno.rds")
temps.in<-readRDS("./GAM Models/ap_egg_temps_in_pheno.rds") 
best.index.phe<-readRDS("./GAM Models/ap_egg_best_index_phe.rds")
thr.geo<-readRDS("./GAM Models/ap_egg_thr_geo.rds")
best.index.geo<-readRDS("./GAM Models/ap_egg_best_index_geo.rds")
vc.pheno<-readRDS("./GAM Models/ap_egg_vc_pheno.rds")
vc.geo<-readRDS("./GAM Models/ap_egg_vc_geo.rds")

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
#generate prediction grid: 
nlat=120
nlon=120
latd=seq(min(apsub$lat),max(apsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(apsub$lon),max(apsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          apsub$lat,apsub$lon)
  grid.extent$dist[k]<-min(dist)
}


# Plot Base Geography Model -----------------------------------------------
grid.extent$year<-as.numeric(2013) #for base, just pick a year that has coverage
grid.extent$doy<-as.numeric(median(apsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(apsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$reg.SST<-mean(apsub$reg.SST,na.rm=TRUE) 
grid.extent$pred<-predict(eg.base,newdata=grid.extent)

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(apsub$lon),ylim=range(apsub$lat),main='Alaska Plaice Distribution, Eggs',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),legend.shrink=0.3)
points(apsub$lon[apsub$Cper10m2==0],apsub$lat[apsub$Cper10m2==0],pch='+',col='white')
symbols(apsub$lon[apsub$Cper10m2>0],
        apsub$lat[apsub$Cper10m2>0],
        circles=log(apsub$Cper10m2+1)[apsub$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)


# Plot Base Phenology Model -----------------------------------------------
grid.extent2<-data.frame('lon'=rep(-155,100),
                         'lat'=rep(51,100),'doy'=seq(min(apsub$doy),max(apsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(apsub$bottom_depth,na.rm=TRUE),100))#setting up another clean grid on which to predict
grid.extent2$pred<-predict(eg.base,newdata=grid.extent2)
grid.extent2$se<-predict(eg.base,newdata=grid.extent2,se=T)[[2]] #select the standard error 
grid.extent2$pred.u<-grid.extent2$pred+1.96*grid.extent2$se #calculate upper confidence interval (95% CI)
grid.extent2$pred.l<-grid.extent2$pred-1.96*grid.extent2$se #calculate lower CI 


windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
plot(grid.extent2$doy,grid.extent2$pred,main='Base Phenology',type='l',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),xlab='Day of Year',cex.lab=1,
     cex.axis=1,cex.main=1,cex.axis=0.9,xlim=range(apsub$doy),
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred.l)),
     col='honeydew2',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col='honeydew2',lty=0)
lines(grid.extent2$doy,grid.extent2$pred,col='grey43')
legend('topleft',legend=c(expression(paste("(log(C/(10m"^2,')+1)')),'95% CI'),
       col=c('grey43','honeydew2'),pch=c(NA,15),lty=c(1,NA),lwd=2,cex=1)
abline(h=0,col='lightslategray',lty=2,lwd=2)

#TEMP EFFECT: Calculate Differences Due to Different Temperature Regimes Based on Best Model --------
#start with threshold geography model to find differences between two predictions to calculate local slopes 
#for Alaska plaice, threshold is 2.064 *C. 
nlat=120
nlon=120
latd=seq(min(apsub$lat),max(apsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(apsub$lon),max(apsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#calculate distance to each positive observation
grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          apsub$lat,apsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-2013
grid.extent$doy<-median(apsub$doy)
grid.extent$reg.SST<-mean(apsub$reg.SST[apsub$reg.SST<2.064]) #threshold temp chosen by AIC values
grid.extent$th<-"TRUE"
grid.extent$bottom_depth<-median(apsub$bottom_depth,na.rm=T)
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)
grid.extent$se<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se #95% CI here
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA #remove predictions that are too far from positive data values
grid.extent$reg.SST<-mean(apsub$reg.SST[apsub$reg.SST>2.064])
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
           xlim=range(apsub$lon),ylim=range(apsub$lat),main='Change in AP Distribution Due to Temperature',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=2,
           legend.lab=expression(paste("(log(C/(10m"^2,')+1)')),
           legend.shrink=0.3) #would prefer to have legend within plot margins, and for all font to be times, but not sure how to do that. 
map("worldHires",fill=T,col="seashell2",add=T)

#now add in the Phenology effect from this model (again, using this model because it produced most deviance explained and lowest AIC): 
grid.extent2<-data.frame('lon'=rep(-155,100),'lat'=rep(51,100),'doy'=seq(min(apsub$doy),max(apsub$doy),length=100),
                         'year'=rep(2013,100),'bottom_depth'=rep(median(apsub$bottom_depth,na.rm=TRUE),100),
                         'reg.SST'=mean(apsub$reg.SST[apsub$reg.SST<2.064]),'th'="TRUE")
grid.extent2$pred<-predict(thr.geo,newdata=grid.extent2)
grid.extent2$se<-predict(thr.geo,newdata=grid.extent2,se=T)[[2]] #select the standard error value from predictions 
grid.extent2$pred.u<-grid.extent2$pred+1.645*grid.extent2$se
grid.extent2$pred.l<-grid.extent2$pred-1.645*grid.extent2$se
grid.extent2$reg.SST<-mean(apsub$reg.SST[apsub$reg.SST>2.064])
grid.extent2$th<-"FALSE"
grid.extent2$pred2<-predict(thr.geo,newdata=grid.extent2)
grid.extent2$se2<-predict(thr.geo,newdata=grid.extent2,se=T)[[2]]
grid.extent2$pred2.u<-grid.extent2$pred2+1.645*grid.extent2$se2
grid.extent2$pred2.l<-grid.extent2$pred2-1.645*grid.extent2$se2

warmcol<-adjustcolor('orangered3',alpha.f=0.3)
coolcol<-adjustcolor('honeydew2',alpha.f=0.8)

plot(grid.extent2$doy,grid.extent2$pred,main='Change in Phenology Due to Two Threshold Conditions',type='l',
     ylim=range(c(grid.extent2$pred.u,grid.extent2$pred2.u,grid.extent2$pred.l,grid.extent2$pred2.l)),
     xlim=range(apsub$doy),col='black',lwd=2,xlab='Day of the Year',
     ylab=expression(paste("(log(C/(10m"^2,')+1)')),cex.lab=1,cex.axis=0.9,cex.main=1)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred.l,rev(grid.extent2$pred.u)),
        col=coolcol,lty=0)
lines(grid.extent2$doy,grid.extent2$pred,col='grey43',lwd=2)
lines(grid.extent2$doy,grid.extent2$pred2,col='indianred3',lwd=2)
polygon(c(grid.extent2$doy,rev(grid.extent2$doy)),c(grid.extent2$pred2.l,rev(grid.extent2$pred2.u)),
        col=warmcol,lty=0)
legend('topleft',legend=c(expression(paste(mu, '(<2.064'^0,' C)')),
                          expression(paste(mu,'(>2.064'^0,' C)'))),
       col=c(coolcol,warmcol),pch=c(15,15),lwd=c(2,2),cex=0.8) #legend could be better

