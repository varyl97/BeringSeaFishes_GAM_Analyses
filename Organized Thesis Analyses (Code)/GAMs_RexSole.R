#Generalized Additive Models: Rex sole
#the following code creates generalized additive models for eggs and larvae of rex sole. 
#these analyses form the basis of my MS thesis. 
#egg data uses an averaged sea surface temperature for the month of March in the Southeastern Bering Sea
#March index was chosen because it is two months before the peak of rex sole CPUE (in May), and thus March 
#conditions are likely more relevant to spawning behavior than temperatures in later months. 
#load egg and larval data: 
rxsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_RxEggs.csv',header=TRUE,check.names=TRUE) #for egg GAMs

rxlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_RxLarv_wCTD.csv',header=TRUE,check.names=TRUE) #for larval GAMs

world<-map_data("world")
BSmap<-world[world$long<(-155)&world$lat>50&world$lat<70,]
#syntax: +geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))

###EGGS: Spawning Behavior 
#the following code generates models to examine flexible 
##Load in local and regional temperature index for March (2 mos before peak egg catch in May) 
reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(rxsub)){
  rxsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==rxsub$year[i]]}

#Base Egg Model: simplest formulation 
eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=rxsub,family=tw(link='log'),method='REML')

summary(eg.base)

windows(width=12,height=8)
plot(eg.base,page=1,scheme=2,
     seWithMean=TRUE,scale=0)

windows(width=12,height=8)
par(mfrow=c(2,2))
vis.gam(eg.base,view=c(c("lon","lat")),too.far=0.025,
        color="heat",plot.type="contour",n.grid=20,
        main='Geographic Linear Predictor (vis.gam), Base Egg')
map("world",fill=T,col="snow4",add=T)
plot(eg.base,select=2,
     xlab='Day of Year',ylab='Estimated Effect, Log(Cper10m2+1)',
     main='Phenologic Smooth')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(eg.base,select=3,shade.col="skyblue3")

windows()
par(mfrow=c(2,2),oma=c(1.5,2,1,1))
gam.check(eg.base)
mtext('gam.check Output for Base Egg Rex CPUE',outer=TRUE,cex=1,line=-0.5)

#Phenological Model with Temperature Threshold 
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]
aic.pheno<-NA*(temps.in)
thr.pheno<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  rxsub$th<-factor(rxsub$reg.SST<=temps.in[i])
  thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+
                        s(lon,lat)+
                        s(bottom_depth,k=5)+
                        s(doy,by=th),
                      data=rxsub,family=tw(link='log'),method='REML')
  aic.pheno[i]<-AIC(thr.pheno[[i]])
}

best.index.phe<-order(aic.pheno)[1]
thr.pheno<-thr.pheno[[best.index.phe]]

windows()
plot(temps.in,aic.pheno,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.pheno)),
     main='Temperature Threshold Flexible Phenology',ylab='AIC Index',
     xlab='Temperature (degC)')
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.phe],lty=2,lwd=2,col='steelblue3')

summary(thr.pheno)

windows(width=12,height=8)
plot(thr.pheno,page=1,
     seWithMean=TRUE,scale=0)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(thr.pheno,select=4,main=paste('Below',round(temps.in[best.index.phe],digits=3),sep=" "),
     seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(thr.pheno,select=3,main=paste('Above',round(temps.in[best.index.phe],digits=3),sep=" "),
     seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)


#Geographical Temperature Threshold Model: 
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]

aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  rxsub$th<-factor(rxsub$reg.SST<=temps.in[i])
  thr.geo[[i]]<-gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                      s(lon,lat,by=th),data=rxsub,
                    family=tw(link='log'),method='REML')
  aic.geo[i]<-AIC(thr.geo[[i]])
}

best.index.geo<-order(aic.geo)[1]
thr.geo<-thr.geo[[best.index.geo]]

windows()
plot(temps.in,aic.geo,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.geo)),
     main='Temperature Threshold Flex Geography',xlab="Temperature (degC)")
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.geo],lty=2,lwd=2,col='steelblue3')

summary(thr.geo)

windows(width=12,height=8)
plot(thr.geo,page=1,scale=0,
     seWithMean=TRUE)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(thr.geo,select=4,scheme=2,too.far=0.025,
     main=paste('Below',round(temps.in[best.index.geo],digits=3),sep=" "),
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)
plot(thr.geo,select=3,scheme=2,too.far=0.025,
     main=paste('Above',round(temps.in[best.index.geo],digits=3),sep=" "),
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)

#Now for Variable-Coefficient Models, beginning with VC Phenology: 
#The variable-coefficient formulation is slightly different from the threshold; 
      #the CPUE is predicted by the inclusion of a covariate that varies in relation to regional temperature,
      #rather than one that varies about a threshold regional temperature value. 
vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                s(doy,by=reg.SST),data=rxsub,family=tw(link='log'),
              method='REML')
summary(vc.pheno)

windows(width=12,height=8)
plot(vc.pheno,
     page=1,scale=0,main='V-C Temp Flexible Phenology, RX Eggs',
     seWithMean=TRUE)

windows(width=16,height=8)
par(mfrow=c(1,2))
plot(vc.pheno,select=2,
     main='V-C Regional Temp, Flexible Phenology',xlab='Day of Year',
     seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(vc.pheno,select=4,
     main='V-C Reg. Temp, Deviation from Avg. Pheno Variation',
     xlab='Day of Year',seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)

#Variable-Coefficient Geography Model: 
vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
              s(lon,lat,by=reg.SST),data=rxsub,family=tw(link='log'),
            method='REML')
summary(vc.geo)

windows(width=12,height=8)
plot(vc.geo,
     page=1,seWithMean=TRUE,main='V-C Flex Geo')

windows(width=15,height=8)
par(mfrow=c(1,2))
plot(vc.geo,select=1,scheme=2,too.far=0.025,
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude',
     main='V-C Rx Egg Flex Geo, Avg. Variation')
map("world",fill=T,col="snow4",add=T)
plot(vc.geo,select=4,scheme=2,too.far=0.025,
     xlab='Longitude',ylab='Latitude',seWithMean=TRUE,
     main='V-C Flex Geo, Deviation from Avg. Variation')
map("world",fill=T,col="snow4",add=T)

##Always save the models to make it much easier to work in the future: 
saveRDS(eg.base,file="./GAM Models/rx_egg_base.rds")
saveRDS(thr.pheno,file="./GAM Models/rx_egg_thr_pheno.rds")
saveRDS(temps.in,file="./GAM Models/rx_egg_temps_in.rds")
saveRDS(aic.pheno,file="./GAM Models/rx_egg_aic_pheno.rds")
saveRDS(best.index.phe,file="./GAM Models/rx_egg_best_index_phe.rds")
saveRDS(thr.geo,file="./GAM Models/rx_egg_thr_geo.rds")
saveRDS(aic.geo,file="./GAM Models/rx_egg_aic_geo.rds")
saveRDS(best.index.geo,file="./GAM Models/rx_egg_best_index_geo.rds")
saveRDS(vc.pheno,file="./GAM Models/rx_egg_vc_pheno.rds")
saveRDS(vc.geo,file="./GAM Models/rx_egg_vc_geo.rds")

eg.base<-readRDS("./GAM Models/rx_egg_base.rds")
thr.pheno<-readRDS("./GAM Models/rx_egg_thr_pheno.rds")
temps.in<-readRDS("./GAM Models/rx_egg_temps_in.rds")
best.index.phe<-readRDS("./GAM Models/rx_egg_best_index_phe.rds")
thr.geo<-readRDS("./GAM Models/rx_egg_thr_geo.rds")
best.index.geo<-readRDS("./GAM Models/rx_egg_best_index_geo.rds")
vc.pheno<-readRDS("./GAM Models/rx_egg_vc_pheno.rds")
vc.geo<-readRDS("./GAM Models/rx_egg_vc_geo.rds")

#Now we can validate best model choice with Akaike Information Criterion: 
aic.base<-AIC(eg.base)
aic.thrph<-AIC(thr.pheno)
aic.thrge<-AIC(thr.geo)
aic.vcph<-AIC(vc.pheno)
aic.vcgeo<-AIC(vc.geo)

aic.rxegg<-data.frame('model'=c('Base','Threshold Pheno','Threshold Geo',
                                'VC Pheno','VC Geo'),
                      'AIC_value'=c(aic.base,aic.thrph,aic.thrge,
                                    aic.vcph,aic.vcgeo))

windows(width=15,height=10)
plot(c(1:5),aic.rxegg$AIC_value,main='AIC Results for Rx Egg Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.rxegg$AIC_value,labels=round(aic.rxegg$AIC_value),pos=c(4,3,3,3,2))
legend("bottomright",legend=c('Base','Threshold Pheno','Threshold Geo',
                             'VC Pheno','VC Geo'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)


###LARVAE: Water Mass Associations, Biogeography
##Now we create larval models, which are all additive but one has a two-dimensional smooth
#these models utilize CTD data
#begin with the larval base model: 
lv.base<-gam((Cper10m2+1)~factor(year)+s(doy)+s(lon,lat)+
               s(bottom_depth,k=5),
             data=rxlarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

windows(width=12,height=8)
par(mfrow=c(2,2))
plot(lv.base,select=1,
     seWithMean=TRUE,scale=0,main='Base Larval Presence GAM, W/O Residuals')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.base,select=2,scheme=2,too.far=0.025,
     
     seWithMean=TRUE,scale=0)
map("world",fill=T,col="snow4",add=T)
plot(lv.base,select=3,
     seWithMean=TRUE,scale=0)
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add in salinity in an additive interaction
lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                  s(bottom_depth,k=5)+
                  s(salinity),data=rxlarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

windows(width=12,height=8)
plot(lv.add.sal,page=1,scale=0,
     seWithMean=TRUE)

windows(width=12,height=8)
plot(lv.add.sal,select=4,
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Salinity',
     xlab='Salinity (PSU)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add in temperature in an additive interaction
lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature),data=rxlarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

windows()
plot(lv.add.temp,page=1,
     seWithMean=TRUE,main='Larval Log(Cper10m2+1) w Temp',scale=0)
    #note: scale=0 allows each individual smooth to have its own y-axis scale

windows(width=12,height=8)
plot(lv.add.temp,select=4,
     seWithMean=TRUE,main='Effect of Temperature',
     xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add both temperature and salinity in individual smooths each 
lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature)+s(salinity),data=rxlarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

windows()
plot(lv.temp.sal,page=1,
     main='Larval Log Presence Temp and Sal',
     seWithMean=TRUE,scale=0)

windows()
par(mfrow=c(1,2))
plot(lv.temp.sal,select=1,seWithMean=TRUE,
     main='Temp and Sal Model')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
par(mfrow=c(1,2))
plot(lv.temp.sal,select=4,
     seWithMean=TRUE,main='Effect of Temp',xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=5,
     seWithMean=TRUE,main='Effect of Salinity',xlab='Salinity (psu)')
abline(h=0,col='sienna3',lty=2,lwd=2)

##Now try temperature and salinity in a two-dimensional, mutually interacting smooth: 
lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
             s(salinity,temperature),data=rxlarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)

windows()
par(mfrow=c(1,2))
plot(lv.2d,select=2,seWithMean=TRUE,
     main='Seasonal Presence, 2D Temp+Sal Model')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.2d,select=1,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
plot(lv.2d,select=4,scheme=TRUE,main='Larval Log Presence, 2D Temp and Sal Effect',
     too.far=0.025,
     xlab='Temperature (degC)',ylab='Salinity (psu)')


#checking based on AIC: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.tempsal<-AIC(lv.temp.sal)
aic.2d<-AIC(lv.2d)

aic.rxlarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Add Sal, Temp','Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.tempsal,aic.2d))

windows()
plot(c(1:5),aic.rxlarv$AIC_value,main='AIC Results for Rx Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.rxlarv$AIC_value,labels=round(aic.rxlarv$AIC_value),pos=c(4,1,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Add Sal, Temp','Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)


# Saved Larval Models -----------------------------------------------------
saveRDS(lv.base,"./GAM Models/rx_larvae_base.rds")
saveRDS(lv.add.sal,"./GAM Models/rx_larvae_addsal.rds")
saveRDS(lv.add.temp,"./GAM Models/rx_larvae_addtemp.rds")
saveRDS(lv.temp.sal,"./GAM Models/rx_larvae_addtempsal.rds")
saveRDS(lv.2d,"./GAM Models/rx_larvae_2d.rds")

lv.base<-readRDS("./GAM Models/rx_larvae_base.rds")
lv.add.sal<-readRDS("./GAM Models/rx_larvae_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/rx_larvae_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/rx_larvae_addtempsal.rds")
lv.2d<-readRDS("./GAM Models/rx_larvae_2d.rds")

