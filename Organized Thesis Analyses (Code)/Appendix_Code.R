
####### Code for Cleaning Raw Data and Creating GAM Formulations: 
#Alaska plaice is used as an example, however these protocols were repeated 
#   for all species. Rex sole had only egg GAM formulations applied, while 
#   northern rock sole, Pacific cod, and yellowfin sole had only larval GAM
#   formulations applied.     

#I. Bring in CTD Data
#First load the CTD data for use in larval dataset later: 
allctd<-read.csv(file="./Environmental Data/All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get an easy year reference

allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,]

#II. Clean Raw Ichthyoplankton Data
#Now we can clean the raw data: 
apeggraw<-read.csv(file='BS_PlaiceEggCatch.csv',header=TRUE,check.names=TRUE)
aplarvraw<-read.csv(file='BS_PlaiceLarvaeCatch.csv',header=TRUE,check.names=TRUE)

apeggraw<-apeggraw[apeggraw$HAUL_ID!='1SS02 81 1 60BON 2',] # remove these observations
aplarvraw<-aplarvraw[aplarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

apeggraw$doy<-as.numeric(mdy.date(apeggraw$MONTH_,apeggraw$DAY_,1960)) #add day of year column 
aplarvraw$doy<-as.numeric(mdy.date(aplarvraw$MONTH_,aplarvraw$DAY_,1960))
apeggraw$id<-paste(apeggraw$CRUISE,apeggraw$LAT,apeggraw$LON,apeggraw$GMT_DATE_TIME,apeggraw$MESH,sep="_") #create unique ID for samples
aplarvraw$id<-paste(aplarvraw$CRUISE,aplarvraw$LAT,aplarvraw$LON,aplarvraw$GMT_DATE_TIME,aplarvraw$MESH,sep="_")

apeggraw[apeggraw$NET==3,] #following lines remove duplicate net types from eggs first
apegg_1_2<-apeggraw[apeggraw$NET<3,]
table(apegg_1_2$NET)
tmp<-table(apegg_1_2$id)
id_2<-names(tmp)[tmp==2]
apegg_1<-apegg_1_2[apegg_1_2$NET==1,]
apegg_2<-apegg_1_2[apegg_1_2$id%in%id_2&apegg_1_2$NET==1,]

apegg_1net<-apegg_1_2[!apegg_1_2$id%in%id_2,]
apegg<-rbind(apegg_1net,apegg_2)
dim(apegg)

#now we can clean net duplicates for larvae
aplarvraw[aplarvraw$NET==3,]
aplarv_1_2<-aplarvraw[aplarvraw$NET<3,]
tmp<-table(aplarv_1_2$id)
id_2<-names(tmp)[tmp==2]
aplarv_2<-aplarv_1_2[aplarv_1_2$id%in%id_2&aplarv_1_2$NET==1,]

aplarv_1net<-aplarv_1_2[!aplarv_1_2$id%in%id_2,]
aplarvae<-rbind(aplarv_1net,aplarv_2)
dim(aplarvae)




#remove stations where primary net =/= to Y
apegg<-apegg[apegg$PRIMARY_NET=='Y',]
aplarvae<-aplarvae[aplarvae$PRIMARY_NET=='Y',]

#add distance to the closest positive catch
apegg$dist<-NA;aplarvae$dist<-NA
for(i in 1:nrow(apegg)){
  apegg$dist[i]<-min(distance.function(apegg$LAT[i],apegg$LON[i],apegg$LAT[apegg$LARVALCATCHPER10M2>0],apegg$LON[apegg$LARVALCATCHPER10M2>0]))
  aplarvae$dist[i]<-min(distance.function(aplarvae$LAT[i], aplarvae$LON[i], aplarvae$LAT[aplarvae$LARVALCATCHPER10M2>0], aplarvae$LON[aplarvae$LARVALCATCHPER10M2>0]))
}

attach(aplarvae)
aplarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(aplarvae)

tmp<-table(aplarvae$YEAR_)
index<-match(aplarvae$YEAR_,names(tmp))
aplarvae$SS<-as.numeric(tmp[index]) #create column with total count of stations sampled per year

attach(apegg)
apegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(apegg)

tmp<-table(apegg$YEAR_)
index<-match(apegg$YEAR_,names(tmp))
apegg$SS<-as.numeric(tmp[index])

apegg$DATE<-NULL
apegg$DATE<-paste(apegg$MONTH_,apegg$DAY_,apegg$YEAR_,sep="/") #create a nice date column

aplarvae$DATE<-NULL
aplarvae$DATE<-paste(aplarvae$MONTH_,aplarvae$DAY_,aplarvae$YEAR_,sep="/")

apsub<-apegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')] #select columns of interest from larger data
apsub<-subset(apsub,doy>99&doy<182) #subset based on histogram of catches across days of year
aplarv<-aplarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]
aplarv<-subset(aplarv,LAT<62)
names(apsub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(aplarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

apsub$SSB<-NA
aplarv$SSB<-NA

SSBdf<-read.csv(file='./Ichthyo Data/SSB_allsp.csv',header=TRUE,check.names=TRUE) #add spawning stock biomass for eggs & larvae
SSBdf$SSB<-as.numeric(SSBdf$SSB)
apSSB<-SSBdf[SSBdf$Species=='alaska plaice',]
apSSB<-apSSB[c('Year','SSB')]

for(i in 1:nrow(aplarv)){
  aplarv$SSB[i]<-apSSB$SSB[apSSB$Year==aplarv$year[i]]
}




for(i in 1:nrow(apsub)){
  apsub$SSB[i]<-apSSB$SSB[apSSB$Year==apsub$year[i]]
}

#now we can link the CTD data to the larval data
aplarv.ctd<-aplarv

allctd$Link_ID<-NA
aplarv.ctd$Link_ID<-NA

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
aplarv.ctd$Link_ID<-paste(aplarv.ctd$CRUISE,aplarv.ctd$STATION,aplarv.ctd$HAUL,sep="_")

aplarv.ctd$temperature<-NA
aplarv.ctd$salinity<-NA
aplarv.ctd$CTD_date<-NA
aplarv.ctd$CTD_link<-NA
aplarv.ctd$CTD_time<-NA

for(i in 1:nrow(aplarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==aplarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    aplarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    aplarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    aplarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    aplarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    aplarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(aplarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
aplarv.ctd<-aplarv.ctd[!aplarv.ctd$temperature=="NaN",]#dim: 2501 24

#make sure that the CTD observations are at similar times as larval observations
aplarv.ctd$date_diff<-NA
aplarv.ctd$date<-parse_date_time(aplarv.ctd$date,orders="mdy") 
aplarv.ctd$CTD_date<-parse_date_time(aplarv.ctd$CTD_date,orders="mdy")
aplarv.ctd$date_diff<-difftime(aplarv.ctd$date,aplarv.ctd$CTD_date,units="hour")
aplarv.ctd<-aplarv.ctd[which(aplarv.ctd$date_diff>(-6)&aplarv.ctd$date_diff<6),]
dim(aplarv.ctd) #2429

#Additional trimming based on biological characteristics: 
#Trim eggs to a depth less than 150m (spawn on the middle shelf between 50-100m)
#Trim larvae to a depth less than 150m (larvae transported to coastal areas after hatching)
apsub<-subset(apsub,bottom_depth<151)
aplarv<-subset(aplarv,bottom_depth<151)
aplarv.ctd<-subset(aplarv.ctd,bottom_depth<151)
aplarv.ctd<-subset(aplarv.ctd,lat<62)

#add regional sea surface temperature (SST) index to ap eggs (apsub) data: 
for(i in 1:nrow(apsub)){
  apsub$reg.SST[i]<-reg.sst$SST[reg.sst$year==apsub$year[i]]
}






#III. GAM Formulations & Some Visualizations - Eggs
#First load in regional SST indices: 
library(mgcv)

reg.sst<-read.csv('./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',
                  header=TRUE,check.names=TRUE)

#First the base egg formulation
eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=apsub,family=tw(link='log'),method='REML')

plot(eg.base,shade=FALSE,page=1,seWithMean=TRUE,scheme=2,scale=0)

summary(eg.base)
AIC(eg.base)

#Then variable-coefficient geography formulation 
vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
              s(lon,lat,by=reg.SST),data=apsub,family=tw(link='log'),
            method='REML')

par(mfrow=c(1,2))
plot(vc.geo,select=1,scheme=2,too.far=0.025,shade=FALSE,
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude',
     main='V-C ap Egg Flex Geo, Avg. Variation')
map("world",fill=T,col="snow4",add=T)
plot(vc.geo,select=4,scheme=2,too.far=0.025,shade=FALSE,
     xlab='Longitude',ylab='Latitude',seWithMean=TRUE,
     main='V-C Flex Geo, Deviation from Avg. Variation')
map("world",fill=T,col="snow4",add=T)

summary(vc.geo)
AIC(vc.geo)

#Variable-coefficient phenology formulation 
vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                s(doy,by=reg.SST),data=apsub,family=tw(link='log'),
              method='REML')

par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5))
plot(vc.pheno,select=2,main='Alaska Plaice VC Phenology, Eggs',seWithMean=TRUE,
     ylim=c(-25,11))
abline(h=0,col='mistyrose4',lty=2,lwd=1.3)
par(oma=c(1,1,1,0.5),mar=c(3,3,3,1.5),new=TRUE)
plot(vc.pheno,select=4,seWithMean=TRUE,shade=TRUE,shade.col=col,ylim=c(-25,11))
legend('topright',legend=c('Flexible Phenology Smooth','Deviation from Avg.Phenology'),
       col=c(NA,col),lwd=c(2,2),cex=0.8)
mtext(c("Day of Year","Anomalies in log(CPUE+1)"),side=c(1,2),line=2.5)

summary(vc.pheno)
AIC(vc.pheno)

#Now for the threshold models; first, threshold phenology
#first create a vector of unique SST indices to test empirically 
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]





#then the for loop to test each SST in a model and ultimately select the best model
aic.pheno<-NA*(temps.in)
thr.pheno<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  apsub$th<-factor(apsub$reg.SST<=temps.in[i])
  thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+
                        s(lon,lat)+
                        s(bottom_depth,k=5)+
                        s(doy,by=th),
                      data=apsub,family=tw(link='log'),method='REML')
  aic.pheno[i]<-AIC(thr.pheno[[i]])
}

best.index.phe<-order(aic.pheno)[1] 
thr.pheno<-thr.pheno[[best.index.phe]] #the model corresponding to the lowest AIC score
temps.in[[best.index.phe]] #identify temperature that served as threshold for best model

summary(thr.pheno)
AIC(thr.pheno)

#now the for loop for the threshold geography model; the temps.in vector is used again
aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
  apsub$th<-factor(apsub$reg.SST<=temps.in[i])
  
  thr.geo[[i]]<-  gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                        s(lon,lat,by=th),data=apsub,
                      family=tw(link='log'),method='REML')
  
  aic.geo[i]<-AIC(thr.geo[[i]])
}

best.index.geo<-order(aic.geo)[1]
thr.geo<-thr.geo[[best.index.geo]]
best.index.geo[[1]]

summary(thr.geo)
AIC(thr.geo)





#The threshold geography is the best model of all five created (based on AIC), so we'll move on 
#   to visualize outputs from this model, starting with geographic distribution above and below 
#   the threshold. 
nlat=120 #create a prediction grid over which we can estimate egg catches (log(cper10m2+1))
nlon=120
latd=seq(min(apsub$lat),max(apsub$lat),length.out=nlat) #center grid over study region 
lond=seq(min(apsub$lon),max(apsub$lon),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          apsub$lat,apsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2007)
grid.extent$doy<-as.numeric(median(apsub$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-median(apsub$bottom_depth,na.rm=TRUE)
grid.extent$reg.SST<-NA
grid.extent$th<-"TRUE"
grid.extent$reg.SST<-mean(apsub$reg.SST,na.rm=TRUE) 
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 
grid.extent$th<-"FALSE"
grid.extent$pred2<-predict(thr.geo,newdata=grid.extent)
grid.extent$pred2[grid.extent$dist>30000]<-NA

#plot grid with catch estimations
par(mai=c(1,1,0.5,0.9),mfrow=c(1,2))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"Lajolla"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=c(-180,-155),ylim=c(52,63),main='',zlim=c(-1.5,8.5),
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-1.8,
           legend.lab=expression(paste("log(C/(10m"^2,')+1)')),legend.shrink=0.4)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T)

image.plot(lond,latd,t(matrix(grid.extent$pred2,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"Lajolla"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=c(-180,-155),ylim=c(52,63),main='',zlim=c(-1.5,8.5),
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-1.8,
           legend.lab=expression(paste("log(C/(10m"^2,')+1)')),legend.shrink=0.4)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T)

#We can do a similar operation to visualize change in egg estimations moving from below to above the threshold
grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')




grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          apsub$lat,apsub$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-2007
grid.extent$doy<-median(apsub$doy)
grid.extent$reg.SST<-mean(apsub$reg.SST[apsub$reg.SST<temps.in[best.index.geo]]) #threshold temp chosen by AIC comparisons
grid.extent$th<-"TRUE"
grid.extent$bottom_depth<-median(apsub$bottom_depth,na.rm=T)
grid.extent$pred<-predict(thr.geo,newdata=grid.extent)
grid.extent$se<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
grid.extent$pred.u<-grid.extent$pred+1.96*grid.extent$se #95% CI here
grid.extent$pred.l<-grid.extent$pred-1.96*grid.extent$se
grid.extent$pred[grid.extent$dist>30000]<-NA #remove predictions that are too far from positive data values
grid.extent$reg.SST<-mean(apsub$reg.SST[apsub$reg.SST>temps.in[best.index.geo]])
grid.extent$th<-"FALSE"
grid.extent$pred2<-predict(thr.geo,newdata=grid.extent)
grid.extent$se2<-predict(thr.geo,newdata=grid.extent,se=T)[[2]]
grid.extent$pred2.u<-grid.extent$pred2+1.96*grid.extent$se
grid.extent$pred2.l<-grid.extent$pred2-1.96*grid.extent$se
grid.extent$diff<-grid.extent$pred2-grid.extent$pred #calculate difference between two regimes

grid.extent$sig.pos<-c(grid.extent$pred2.l>grid.extent$pred.u) #isolate areas where there is a higher predicted CPUE at a higher temperature
grid.extent$sig.neg<-c(grid.extent$pred2.u<grid.extent$pred.l)
grid.extent$pos.diff<-grid.extent$diff*grid.extent$sig.pos 
grid.extent$neg.diff<-grid.extent$diff*grid.extent$sig.neg
max.slope<-max(grid.extent$diff,na.rm=T)

#and plot outcome: 
par(mai=c(1,1,0.5,0.5))
image.plot(lond,latd,t(matrix(grid.extent$diff,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=hcl.colors(100,"BuPu",rev=T),ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')), 
           xlim=c(-180,-155),ylim=c(52,63),
           main="",
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("Change in log(C/(10m"^2,')+1)')),
           legend.shrink=0.4)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T)





#IV. GAM Formulations & Some Visualizations - Larvae
# The larval models are mostly additive, with less complicated formulations
lv.base<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
               s(bottom_depth,k=5),
             data=aplarv.ctd,family=tw(link='log'),method='REML')

summary(lv.base)
AIC(lv.base)

#then additive sea surface salinity (SSS)
lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                  s(bottom_depth,k=5)+
                  s(salinity),data=aplarv.ctd,family=tw(link='log'),
                method='REML')

summary(lv.add.sal)
AIC(lv.add.sal)

#then additive SST 
lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature),data=aplarv.ctd,family=tw(link='log'),
                 method='REML')

summary(lv.add.sal)
AIC(lv.add.sal)

#now we add SSS and SST as additive terms
lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature)+s(salinity),data=aplarv.ctd,
                 family=tw(link='log'),method='REML')

summary(lv.add.temp)
AIC(lv.add.temp)

#now we'll create the final model, with a tensor smooth bivariate term for SSS and SST
lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth,k=5)+
             te(salinity,temperature),data=aplarv.ctd,family=tw(link='log'),
           method='REML')

summary(lv.2d)
AIC(lv.2d) #lowest AIC of all models, best model for larvae

#To visualize the best model, the bivariate tensor smooth model, we can use the following code: 
nlat=120 #make a similar prediction grid as that for eggs
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
grid.extent$temperature<-as.numeric(mean(aplarv.ctd$temperature))
grid.extent$salinity<-as.numeric(mean(aplarv.ctd$salinity))
grid.extent$pred<-predict(lv.2d,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA 

par(mai=c(1,1,0.5,0.9),mfrow=c(2,1))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"Lajolla"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=c(-180,-155),ylim=c(52,63),
           main='',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2,
           legend.lab=expression(paste("log(C/(10m"^2,')+1)')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T)





























