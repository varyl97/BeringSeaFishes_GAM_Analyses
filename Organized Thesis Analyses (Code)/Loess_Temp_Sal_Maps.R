## The following document creates loess interpolated salinity and 
      #temperature maps for the Bering Sea across years with substantial 
      #observations. 
## these interpolations utilize CTD data before it is truncated to match 
      #ichthyoplankton data coverage. Temperature and salinity values are for 
      #top 10 meters of the water column. 
allctd<-read.csv(file="./Environmental Data/All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                                 'Temperature','Conductivity','Salinity',
                 'Sigma.T','Flag','Cruise','Station','Haul','Grid',
                 'FOCI_HAUL_ID','FOCI_file')]
 
allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference
allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,]
newctd<-allctd[allctd$Depth<11,]

newctd$ider<-paste0(newctd$Cruise,newctd$Haul,newctd$Station,newctd$Date)
unique<-unique(newctd$ider)

uniqname<-paste0("sys",newctd$ider)

splitCTD<-split(newctd,uniqname)

mean_function<-function(x){
  x$mean_temp<-mean(x$Temperature)
  x$mean_sal<-mean(x$Salinity)
  return(x)}

splitCTD2<-lapply(splitCTD,mean_function)

library(rlist)
allctd2<-list.rbind(splitCTD2)

#Start here: 
muctd<-read.csv(file="./Environmental Data/Split_CTD_MeanTopTempSal.csv",
                header=T,check.names=T)


################# FOR THESIS: #######################################
#average SST and SSS for region: 
## fix this so that it's just the months of interest
#cutting just to may
muctd$month<-str_sub(muctd$Date,start=1,end=1) #not perfect for double digit months but doesn't matter 
cutctd<-muctd[muctd$month==5,]

cutctd<-cutctd[!is.na(cutctd$Longitude),]
cutctd<-cutctd[!is.na(cutctd$Latitude),]

lon.poly<-c(-154.2265, -166.4125, -161.5889, -154.2688)
lat.poly<-c(58.32915, 53.87126, 53.12828, 55.95868)

point.in.poly<-in.chull(cutctd$Longitude,cutctd$Latitude,lon.poly,lat.poly)
cutctd<-cutctd[!point.in.poly==TRUE,]

zlim.t<-range(cutctd$mean_temp,na.rm=T)
zlim.s<-range(cutctd$mean_sal,na.rm=T)

loess.t<-loess(mean_temp~Longitude*Latitude,
               data=cutctd,degree=2,span=0.15)
loess.s<-loess(mean_sal~Longitude*Latitude,
               data=cutctd,degree=2,span=0.15)

lond<-seq(min(cutctd$Longitude,na.rm=T),max(cutctd$Longitude,na.rm=T),
          length.out=100)
latd<-seq(min(cutctd$Latitude,na.rm=T),max(cutctd$Latitude,na.rm=T),
          length.out=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist[j]<-min(distance.function(
    predict.grid$Latitude[j],predict.grid$Longitude[j],
    cutctd$Latitude,cutctd$Longitude),na.rm=T)}
predict.grid[predict.grid$dist>30000,]<-NA

sst.pred<-predict(loess.t,newdata=predict.grid)
sss.pred<-predict(loess.s,newdata=predict.grid)

windows(width=12,height=5);par(mai=c(1,1,0.5,1.2),mfrow=c(1,2))
image.plot(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),
           zlim=zlim.t,ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Temperature ("^0,"C)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond,latd,sst.pred,levels=c(2),add=T,lwd=1.9,col="grey90")
contour(lond,latd,sst.pred,levels=c(6),add=T,lwd=1.9,col='grey90')
map("worldHires",fill=T,col="gainsboro",add=T)

image.plot(lond,latd,sss.pred,col=hcl.colors(100,"Mako",rev=T),
           zlim=zlim.s,ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Salinity (psu)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond,latd,sss.pred,levels=c(31),add=T,lwd=1.9,col="grey90")
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=1.9,col='grey90')
map("worldHires",fill=T,col="gainsboro",add=T)

#averages in warm and cold regimes: 
#cold regime: '06 through '13
#warm regime: '14 through '16

cutctd<-cutctd %>% mutate(regime = 
                            case_when(Year>=2006 & Year<=2013 ~ "Cold",
                                      Year>=2001 & Year<=2005 ~ "Warm",
                                      Year>=2014 & Year<=2016 ~ "Warm"))

warm<-cutctd[cutctd$regime=="Warm",]
cold<-cutctd[cutctd$regime=="Cold",]

zlimt.cold<-range(cold$mean_temp,na.rm=T)
zlimt.warm<-range(warm$mean_temp,na.rm=T)

zlims.cold<-range(cold$mean_sal,na.rm=T)
zlims.warm<-range(warm$mean_sal,na.rm=T)

loess.t.cold<-loess(mean_temp~Longitude*Latitude,
                    data=cold,degree=2,span=0.15)
loess.t.warm<-loess(mean_temp~Longitude*Latitude,
                    data=warm,degree=2,span=0.15)

loess.s.cold<-loess(mean_sal~Longitude*Latitude,
                    data=cold,degree=2,span=0.15)
loess.s.warm<-loess(mean_sal~Longitude*Latitude,
                  data=warm,degree=2,span=0.15)

lond.cold<-seq(min(cold$Longitude,na.rm=T),
               max(cold$Longitude,na.rm=T),length=100)
latd.cold<-seq(min(cold$Latitude,na.rm=T),
               max(cold$Latitude,na.rm=T),length=100)
lond.warm<-seq(min(warm$Longitude,na.rm=T),
               max(warm$Longitude,na.rm=T),length=100)
latd.warm<-seq(min(warm$Latitude,na.rm=T),
               max(warm$Latitude,na.rm=T),length=100)

grid.cold<-expand.grid(lond.cold,latd.cold)
names(grid.cold)<-c("Longitude","Latitude")

grid.warm<-expand.grid(lond.warm,latd.warm)
names(grid.warm)<-c("Longitude","Latitude")

for(j in 1:nrow(grid.cold)){
    grid.cold$dist[j]<-min(distance.function(
    grid.cold$Latitude[j],grid.cold$Longitude[j],
    cold$Latitude,cold$Longitude),na.rm=T)}
grid.cold[grid.cold$dist>30000,]<-NA

for(j in 1:nrow(grid.warm)){
    grid.warm$dist[j]<-min(distance.function(
    grid.warm$Latitude[j],grid.warm$Longitude[j],
    warm$Latitude,warm$Longitude),na.rm=T)}
grid.warm[grid.warm$dist>30000,]<-NA

sst.cold<-predict(loess.t.cold,newdata=grid.cold)
sst.warm<-predict(loess.t.warm,newdata=grid.warm)

sss.cold<-predict(loess.s.cold,newdata=grid.cold)
sss.warm<-predict(loess.s.warm,newdata=grid.warm)

windows(width=12,height=5);par(mai=c(1,1,0.5,1.2),mfrow=c(1,2))
image.plot(lond.warm,latd.warm,sst.warm,col=hcl.colors(100,"Viridis"),
           zlim=c(-2,9.5),ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Temperature ("^0,"C)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond.warm,latd.warm,sst.warm,levels=c(2),add=T,lwd=1.75,col="grey")
contour(lond.warm,latd.warm,sst.warm,levels=c(6),add=T,lwd=1.75,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

image.plot(lond.warm,latd.warm,sss.warm,col=hcl.colors(100,"Mako",rev=T),
           zlim=c(29.5,33.5),ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Salinity (psu)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond.warm,latd.warm,sss.warm,levels=c(31),add=T,lwd=1.75,col="grey")
contour(lond.warm,latd.warm,sss.warm,levels=c(32),add=T,lwd=1.75,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

windows(width=12,height=5);par(mai=c(1,1,0.5,1.2),mfrow=c(1,2))
image.plot(lond.cold,latd.cold,sst.cold,col=hcl.colors(100,"Viridis"),
           zlim=c(-2,9.5),ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Temperature ("^0,"C)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond.cold,latd.cold,sst.cold,levels=c(2),add=T,lwd=1.75,col="grey")
contour(lond.cold,latd.cold,sst.cold,levels=c(6),add=T,lwd=1.75,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

image.plot(lond.cold,latd.cold,sss.cold,col=hcl.colors(100,"Mako",rev=T),
           zlim=c(29.5,33.5),ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           main="",
           xlim=c(-174,-155),ylim=c(53,60),
           legend.lab=expression(paste("Salinity (psu)")),
           legend.line=(-1.8))
contour(bathy,levels=-c(50,100),labcex=0.4,col='grey58',add=T)
contour(lond.cold,latd.cold,sss.cold,levels=c(31),add=T,lwd=1.75,col="grey")
contour(lond.cold,latd.cold,sss.cold,levels=c(32),add=T,lwd=1.75,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)






###by year (some older code below here)
year<-c(2003,2005,2008,2009) #2 cold and warm years, respectively, according to Baker 202102.
#cold 1=2003, 2=2005; warm 1=2008, warm 2=2009

for(i in 1:length(year)){
  sub<-muCTD[muCTD$Year==year[i],]
}

zlim<-range(sub$mean_temp,na.rm=T)
sub<-sub[!is.na(sub$mean_temp),]

cold1<-muCTD[muCTD$Year==2003,]
zlim1<-range(cold1$mean_temp,na.rm=T)
cold2<-muCTD[muCTD$Year==2005,]
zlim2<-range(cold2$mean_temp,na.rm=T)

warm1<-muCTD[muCTD$Year==2008,]
zlim3<-range(warm1$mean_temp,na.rm=T)
warm2<-muCTD[muCTD$Year==2009,]
zlim4<-range(warm2$mean_temp,na.rm=T)



#temperature plots: 
loess.t1<-loess(mean_temp~Longitude*Latitude,
                data=cold1,span=0.15,degree=2)
loess.t2<-loess(mean_temp~Longitude*Latitude,
                data=cold2,span=0.15,degree=2)
loess.t3<-loess(mean_temp~Longitude*Latitude,
                data=warm1,span=0.15,degree=2)
loess.t4<-loess(mean_temp~Longitude*Latitude,
                data=warm2,span=0.15,degree=2)

#framework of prediction grid
lond<-seq(min(muCTD$Longitude,na.rm=T),max(muCTD$Longitude,na.rm=T),
          length=100)
latd<-seq(min(muCTD$Latitude,na.rm=T),max(muCTD$Latitude,na.rm=T),
          length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

lond1<-seq(min(cold1$Longitude,na.rm=T),
          max(cold1$Longitude,na.rm=T),length=100)
latd1<-seq(min(cold1$Latitude,na.rm=T),
          max(cold1$Latitude,na.rm=T),length=100)
predict.grid1<-expand.grid(lond1,latd1)
names(predict.grid1)<-c("Longitude","Latitude")

lond2<-seq(min(cold2$Longitude,na.rm=T),
           max(cold2$Longitude,na.rm=T),length=100)
latd2<-seq(min(cold2$Latitude,na.rm=T),
           max(cold2$Latitude,na.rm=T),length=100)
predict.grid2<-expand.grid(lond2,latd2)
names(predict.grid2)<-c("Longitude","Latitude")

lond3<-seq(min(warm1$Longitude,na.rm=T),
           max(warm1$Longitude,na.rm=T),length=100)
latd3<-seq(min(warm1$Latitude,na.rm=T),
           max(warm1$Latitude,na.rm=T),length=100)
predict.grid3<-expand.grid(lond3,latd3)
names(predict.grid3)<-c("Longitude","Latitude")

lond4<-seq(min(warm2$Longitude,na.rm=T),
           max(warm2$Longitude,na.rm=T),length=100)
latd4<-seq(min(warm2$Latitude,na.rm=T),
           max(warm2$Latitude,na.rm=T),length=100)
predict.grid4<-expand.grid(lond4,latd4)
names(predict.grid4)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist[j]<-min(distance.function(
    predict.grid$Latitude[j],predict.grid$Longitude[j],
    muCTD$Latitude[j],muCTD$Longitude[j]))}

for(j in 1:nrow(predict.grid1)){
  predict.grid1$dist[j]<-min(distance.function(
    predict.grid1$Latitude[j],predict.grid1$Longitude[j],
    cold1$Latitude[j],cold1$Longitude[j]))}

for(j in 1:nrow(predict.grid2)){
  predict.grid2$dist[j]<-min(distance.function(
    predict.grid2$Latitude[j],predict.grid2$Longitude[j],
    cold2$Latitude[j],cold2$Longitude[j]))}

for(j in 1:nrow(predict.grid3)){
  predict.grid3$dist[j]<-min(distance.function(
    predict.grid3$Latitude[j],predict.grid3$Longitude[j],
    warm2$Latitude[j],warm1$Longitude[j]))}

for(j in 1:nrow(predict.grid4)){
  predict.grid4$dist[j]<-min(distance.function(
    predict.grid4$Latitude[j],predict.grid4$Longitude[j],
    warm2$Latitude[j],warm2$Longitude[j]))}

sst.pred<-predict(loess.all,newdata=predict.grid)

sst.pred1<-predict(loess.t1,newdata=predict.grid1)
#sst.pred1[predict.grid1$dist>500000]<-NA

sst.pred2<-predict(loess.t2,newdata=predict.grid2)
#sst.pred2[predict.grid2$dist>500000]<-NA

sst.pred3<-predict(loess.t3,newdata=predict.grid3)
#sst.pred3[predict.grid3$dist>500000]<-NA

sst.pred4<-predict(loess.t4,newdata=predict.grid4)
#sst.pred4[predict.grid4$dist>500000]<-NA

windows(width=6,height=5);par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),
      zlim=zlim,ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),
      main="",
      xlim=c(-175,-155),ylim=c(53,65),
      legend.lab=expression(paste("Temperature ("^0,"C)")),
      legend.line=(-2.2))
contour(lond,latd,sst.pred,levels=c(2),add=T,lwd=1.5,col="grey")
contour(lond,latd,sst.pred,levels=c(6),add=T,lwd=1.5,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

windows(width=15,height=13)
par(mfrow=c(2,2))
image(lond1,latd1,sst.pred1,col=hcl.colors(100,"Lajolla"),
      zlim=zlim,ylab='Latitude',xlab='Longitude',
      main=paste('Top 10m Temp, 2003'),
      xlim=c(-177,-160),ylim=c(54,63))
contour(lond1,latd1,sst.pred1,levels=c(2.5),add=T,lwd=2)
contour(lond1,latd1,sst.pred1,levels=c(-1),add=T,lwd=2,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

image(lond2,latd2,sst.pred2,col=hcl.colors(100,"Lajolla"),
      zlim=zlim,ylab='Latitude',xlab='Longitude',
      main=paste('Top 10m Temp, 2004'),
      xlim=c(-177,-160),ylim=c(54,63))
contour(lond1,latd1,sst.pred1,levels=c(2.5),add=T,lwd=2)
contour(lond1,latd1,sst.pred1,levels=c(-1),add=T,lwd=2,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

image(lond3,latd3,sst.pred1,col=hcl.colors(100,"Lajolla"),
      zlim=zlim,ylab='Latitude',xlab='Longitude',
      main=paste('Top 10m Temp, 2008'),
      xlim=c(-177,-160),ylim=c(54,63))
contour(lond3,latd3,sst.pred3,levels=c(2.5),add=T,lwd=2)
contour(lond3,latd3,sst.pred3,levels=c(-1),add=T,lwd=2,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)

image(lond4,latd4,sst.pred4,col=hcl.colors(100,"Lajolla"),
      zlim=zlim,ylab='Latitude',xlab='Longitude',
      main=paste('Top 10m Temp, 2009'),
      xlim=c(-177,-160),ylim=c(54,63))
contour(lond4,latd4,sst.pred4,levels=c(2.5),add=T,lwd=2)
contour(lond4,latd4,sst.pred4,levels=c(-1),add=T,lwd=2,col='grey')
map("worldHires",fill=T,col="gainsboro",add=T)


###### For thesis: 























#monthly iteration - 2/21/2022
#still really patchy, think the data just might not be ideal 
muctd<-muCTD
head(muctd)

#2006
may<-muctd[muctd$Month==5,]
may06<-may[may$Year==2006,] #had most observations
may06<-may06[!is.na(may06$mean_temp),]

zlim<-range(may06$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=may06,span=0.35,degree=2)
lond<-seq(min(may06$Longitude,na.rm=T),max(may06$Longitude,na.rm=T),length=100)
latd<-seq(min(may06$Latitude,na.rm=T),max(may06$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           may06$Latitude[j],
                                           may06$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=15,height=13)
par(mfrow=c(2,2))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab='Latitude',xlab='Longitude',main='May 2006',
      xlim=c(-175,-158),ylim=c(50,60))
contour(lond,latd,sst.pred,levels=c(-1.5),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
mtext("Sea Surface Mean Temperature",outer=TRUE,line=-1)

#2007
may07<-may[may$Year==2007,] #had most observations
may07<-may07[!is.na(may07$mean_temp),]

zlim<-range(may07$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=may07,span=0.35,degree=2)
lond<-seq(min(may07$Longitude,na.rm=T),max(may07$Longitude,na.rm=T),length=100)
latd<-seq(min(may07$Latitude,na.rm=T),max(may07$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           may07$Latitude[j],
                                           may07$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab='Latitude',xlab='Longitude',main='May 2007',
      xlim=c(-175,-158),ylim=c(50,60))
contour(lond,latd,sst.pred,levels=c(-1.5),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)

#2014
may14<-may[may$Year==2014,] 
may14<-may14[!is.na(may14$mean_temp),]

zlim<-range(may14$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=may14,span=0.35,degree=2)
lond<-seq(min(may14$Longitude,na.rm=T),max(may14$Longitude,na.rm=T),length=100)
latd<-seq(min(may14$Latitude,na.rm=T),max(may14$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           may14$Latitude[j],
                                           may14$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab='Latitude',xlab='Longitude',main='May 2014',
      xlim=c(-175,-158),ylim=c(50,60))
contour(lond,latd,sst.pred,levels=c(-1.5),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)

#2016
may16<-may[may$Year==2016,] 
may16<-may16[!is.na(may16$mean_temp),]

zlim<-range(may16$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=may16,span=0.35,degree=2)
lond<-seq(min(may16$Longitude,na.rm=T),max(may16$Longitude,na.rm=T),length=100)
latd<-seq(min(may16$Latitude,na.rm=T),max(may16$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           may16$Latitude[j],
                                           may16$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab='Latitude',xlab='Longitude',main='May 2016',
      xlim=c(-175,-158),ylim=c(50,60))
contour(lond,latd,sst.pred,levels=c(-1.5),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)

###### Spring & Summer months iteration: Apr - May - June - July: 
#years: 1997, 2000, 2002, 2005-2017 
#sea surface temperature: 

muctd<-read.csv(file="./Environmental Data/Split_CTD_MeanTopTempSal.csv",
                header=T,check.names=T)
head(muctd)
muctd$Month<-str_sub(muctd$Date,start=1,end=1)

#for plotting: 
str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name)
bathybath<-as.bathy(bathy)


sprsum<-muctd[muctd$Month>3&muctd$Month<8,]

#1997
spr97<-sprsum[sprsum$Year==1997,] 
spr97<-spr97[!is.na(spr97$mean_temp),]

zlim<-range(spr97$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr97,span=0.35,degree=2)
lond<-seq(min(spr97$Longitude,na.rm=T),max(spr97$Longitude,na.rm=T),length=100)
latd<-seq(min(spr97$Latitude,na.rm=T),max(spr97$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr97$Latitude[j],
                                           spr97$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='1997',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2000
spr00<-sprsum[sprsum$Year==2000,] 
spr00<-spr00[!is.na(spr00$mean_temp),]

zlim<-range(spr00$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr00,span=0.35,degree=2)
lond<-seq(min(spr00$Longitude,na.rm=T),max(spr00$Longitude,na.rm=T),length=100)
latd<-seq(min(spr00$Latitude,na.rm=T),max(spr00$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr00$Latitude[j],
                                           spr00$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2000',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col="grey73",add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)


#2002
spr02<-sprsum[sprsum$Year==2002,] 
spr02<-spr02[!is.na(spr02$mean_temp),]

zlim<-range(spr02$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr02,span=0.35,degree=2)
lond<-seq(min(spr02$Longitude,na.rm=T),max(spr02$Longitude,na.rm=T),length=100)
latd<-seq(min(spr02$Latitude,na.rm=T),max(spr02$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr02$Latitude[j],
                                           spr02$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2002',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),col="grey73",add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2005
spr05<-sprsum[sprsum$Year==2005,] 
spr05<-spr05[!is.na(spr05$mean_temp),]

zlim<-range(spr05$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr05,span=0.35,degree=2)
lond<-seq(min(spr05$Longitude,na.rm=T),max(spr05$Longitude,na.rm=T),length=100)
latd<-seq(min(spr05$Latitude,na.rm=T),max(spr05$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr05$Latitude[j],
                                           spr05$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2005',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col="grey73",add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2006
spr06<-sprsum[sprsum$Year==2006,] 
spr06<-spr06[!is.na(spr06$mean_temp),]

zlim<-range(spr06$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr06,span=0.35,degree=2)
lond<-seq(min(spr06$Longitude,na.rm=T),max(spr06$Longitude,na.rm=T),length=100)
latd<-seq(min(spr06$Latitude,na.rm=T),max(spr06$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr06$Latitude[j],
                                           spr06$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2006',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2007
spr07<-sprsum[sprsum$Year==2007,] 
spr07<-spr07[!is.na(spr07$mean_temp),]

zlim<-range(spr07$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr07,span=0.35,degree=2)
lond<-seq(min(spr07$Longitude,na.rm=T),max(spr07$Longitude,na.rm=T),length=100)
latd<-seq(min(spr07$Latitude,na.rm=T),max(spr07$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr07$Latitude[j],
                                           spr07$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2007',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col="grey73",add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2008
spr08<-sprsum[sprsum$Year==2008,] 
spr08<-spr08[!is.na(spr08$mean_temp),]

zlim<-range(spr08$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr08,span=0.35,degree=2)
lond<-seq(min(spr08$Longitude,na.rm=T),max(spr08$Longitude,na.rm=T),length=100)
latd<-seq(min(spr08$Latitude,na.rm=T),max(spr08$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr08$Latitude[j],
                                           spr08$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2008',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2009
spr09<-sprsum[sprsum$Year==2009,] #had most observations
spr09<-spr09[!is.na(spr09$mean_temp),]

zlim<-range(spr09$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr09,span=0.35,degree=2)
lond<-seq(min(spr09$Longitude,na.rm=T),max(spr09$Longitude,na.rm=T),length=100)
latd<-seq(min(spr09$Latitude,na.rm=T),max(spr09$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr09$Latitude[j],
                                           spr09$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2009',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2010
spr10<-sprsum[sprsum$Year==2010,] #had most observations
spr10<-spr10[!is.na(spr10$mean_temp),]

zlim<-range(spr10$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr10,span=0.35,degree=2)
lond<-seq(min(spr10$Longitude,na.rm=T),max(spr10$Longitude,na.rm=T),length=100)
latd<-seq(min(spr10$Latitude,na.rm=T),max(spr10$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr10$Latitude[j],
                                           spr10$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2010',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2011
spr11<-sprsum[sprsum$Year==2011,] #had most observations
spr11<-spr11[!is.na(spr11$mean_temp),]

zlim<-range(spr11$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr11,span=0.35,degree=2)
lond<-seq(min(spr11$Longitude,na.rm=T),max(spr11$Longitude,na.rm=T),length=100)
latd<-seq(min(spr11$Latitude,na.rm=T),max(spr11$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr11$Latitude[j],
                                           spr11$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2011',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2012
spr12<-sprsum[sprsum$Year==2012,] #had most observations
spr12<-spr12[!is.na(spr12$mean_temp),]

zlim<-range(spr12$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr12,span=0.35,degree=2)
lond<-seq(min(spr12$Longitude,na.rm=T),max(spr12$Longitude,na.rm=T),length=100)
latd<-seq(min(spr12$Latitude,na.rm=T),max(spr12$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr12$Latitude[j],
                                           spr12$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2012',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2013
spr13<-sprsum[sprsum$Year==2013,] #had most observations
spr13<-spr13[!is.na(spr13$mean_temp),]

zlim<-range(spr13$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr13,span=0.35,degree=2)
lond<-seq(min(spr13$Longitude,na.rm=T),max(spr13$Longitude,na.rm=T),length=100)
latd<-seq(min(spr13$Latitude,na.rm=T),max(spr13$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr13$Latitude[j],
                                           spr13$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2013',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2014
spr14<-sprsum[sprsum$Year==2014,] 
spr14<-spr14[!is.na(spr14$mean_temp),]

zlim<-range(spr14$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr14,span=0.35,degree=2)
lond<-seq(min(spr14$Longitude,na.rm=T),max(spr14$Longitude,na.rm=T),length=100)
latd<-seq(min(spr14$Latitude,na.rm=T),max(spr14$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr14$Latitude[j],
                                           spr14$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2014',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2015
spr15<-sprsum[sprsum$Year==2015,] #had most observations
spr15<-spr15[!is.na(spr15$mean_temp),]

zlim<-range(spr15$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr15,span=0.35,degree=2)
lond<-seq(min(spr15$Longitude,na.rm=T),max(spr15$Longitude,na.rm=T),length=100)
latd<-seq(min(spr15$Latitude,na.rm=T),max(spr15$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr15$Latitude[j],
                                           spr15$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2015',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2016
spr16<-sprsum[sprsum$Year==2016,] 
spr16<-spr16[!is.na(spr16$mean_temp),]

zlim<-range(spr16$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr16,span=0.35,degree=2)
lond<-seq(min(spr16$Longitude,na.rm=T),max(spr16$Longitude,na.rm=T),length=100)
latd<-seq(min(spr16$Latitude,na.rm=T),max(spr16$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr16$Latitude[j],
                                           spr16$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2016',
      xlim=c(-180,-155),ylim=c(52,63))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

#2017
spr17<-sprsum[sprsum$Year==2017,] #had most observations
spr17<-spr17[!is.na(spr17$mean_temp),]

zlim<-range(spr17$mean_temp,na.rm=T)

loess<-loess(mean_temp~Longitude*Latitude,data=spr17,span=0.35,degree=2)
lond<-seq(min(spr17$Longitude,na.rm=T),max(spr17$Longitude,na.rm=T),length=100)
latd<-seq(min(spr17$Latitude,na.rm=T),max(spr17$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr17$Latitude[j],
                                           spr17$Longitude[j]))}

sst.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sst.pred,col=hcl.colors(100,"Viridis"),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2017',
      xlim=c(-180,-155),ylim=c(60,70))
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey73',add=T)
contour(lond,latd,sst.pred,levels=c(0),add=T,lwd=2)
contour(lond,latd,sst.pred,levels=c(3),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"Viridis"),
           legend.shrink=0.3,legend.width=1)
mtext("SST (*C)",side=4,line=0.75)

### Mean Salinity: 

#years: 1997, 2000, 2002, 2005-2017 
muctd<-muCTD
head(muctd)
sprsum<-muctd[muctd$Month>3&muctd$Month<8,]

#1997
spr97<-sprsum[sprsum$Year==1997,] 
spr97<-spr97[!is.na(spr97$mean_sal),]

zlim<-range(spr97$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr97,span=0.35,degree=2)
lond<-seq(min(spr97$Longitude,na.rm=T),max(spr97$Longitude,na.rm=T),length=100)
latd<-seq(min(spr97$Latitude,na.rm=T),max(spr97$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr97$Latitude[j],
                                           spr97$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='1997',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2000
spr00<-sprsum[sprsum$Year==2000,] 
spr00<-spr00[!is.na(spr00$mean_sal),]

zlim<-range(spr00$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr00,span=0.35,degree=2)
lond<-seq(min(spr00$Longitude,na.rm=T),max(spr00$Longitude,na.rm=T),length=100)
latd<-seq(min(spr00$Latitude,na.rm=T),max(spr00$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr00$Latitude[j],
                                           spr00$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2000',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2002
spr02<-sprsum[sprsum$Year==2002,] 
spr02<-spr02[!is.na(spr02$mean_sal),]

zlim<-range(spr02$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr02,span=0.35,degree=2)
lond<-seq(min(spr02$Longitude,na.rm=T),max(spr02$Longitude,na.rm=T),length=100)
latd<-seq(min(spr02$Latitude,na.rm=T),max(spr02$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr02$Latitude[j],
                                           spr02$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2002',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2005
spr05<-sprsum[sprsum$Year==2005,] 
spr05<-spr05[!is.na(spr05$mean_sal),]

zlim<-range(spr05$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr05,span=0.35,degree=2)
lond<-seq(min(spr05$Longitude,na.rm=T),max(spr05$Longitude,na.rm=T),length=100)
latd<-seq(min(spr05$Latitude,na.rm=T),max(spr05$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr05$Latitude[j],
                                           spr05$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2005',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2006
spr06<-sprsum[sprsum$Year==2006,] 
spr06<-spr06[!is.na(spr06$mean_sal),]

zlim<-range(spr06$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr06,span=0.35,degree=2)
lond<-seq(min(spr06$Longitude,na.rm=T),max(spr06$Longitude,na.rm=T),length=100)
latd<-seq(min(spr06$Latitude,na.rm=T),max(spr06$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr06$Latitude[j],
                                           spr06$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=18,height=15)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2006',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2007
spr07<-sprsum[sprsum$Year==2007,] 
spr07<-spr07[!is.na(spr07$mean_sal),]

zlim<-range(spr07$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr07,span=0.35,degree=2)
lond<-seq(min(spr07$Longitude,na.rm=T),max(spr07$Longitude,na.rm=T),length=100)
latd<-seq(min(spr07$Latitude,na.rm=T),max(spr07$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr07$Latitude[j],
                                           spr07$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2007',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2008 - didn't give anything useful... lack of data
#spr08<-sprsum[sprsum$Year==2008,] 
#spr08<-spr08[!is.na(spr08$mean_sal),]

#zlim<-range(spr08$mean_sal,na.rm=T)

#loess<-loess(mean_sal~Longitude*Latitude,data=spr08,span=0.35,degree=2)
#lond<-seq(min(spr08$Longitude,na.rm=T),max(spr08$Longitude,na.rm=T),length=100)
#latd<-seq(min(spr08$Latitude,na.rm=T),max(spr08$Latitude,na.rm=T),length=100)
#predict.grid<-expand.grid(lond,latd)
#names(predict.grid)<-c("Longitude","Latitude")

#for(j in 1:nrow(predict.grid)){
 # predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                #                           predict.grid$Longitude[j],
                #                           spr08$Latitude[j],
                #                           spr08$Longitude[j]))}

#sss.pred<-predict(loess,newdata=predict.grid)

#windows(width=6,height=5)
#par(oma=c(1,1,1,4))
#image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
 #     ylab=expression(paste("Latitude ("^0,'N)')),
  #    xlab=expression(paste("Longitude ("^0,'E)')),main='2008',
   #   xlim=c(-180,-155),ylim=c(52,63))
#contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
#map("worldHires",fill=T,col="gainsboro",add=T)
#par(oma=c(1,1,1,1.5))
#image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
 #          legend.shrink=0.3,legend.width=1)
#mtext("SSS (psu)",side=4,line=1.25)

#2009
spr09<-sprsum[sprsum$Year==2009,] #had most observations
spr09<-spr09[!is.na(spr09$mean_sal),]

zlim<-range(spr09$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr09,span=0.35,degree=2)
lond<-seq(min(spr09$Longitude,na.rm=T),max(spr09$Longitude,na.rm=T),length=100)
latd<-seq(min(spr09$Latitude,na.rm=T),max(spr09$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr09$Latitude[j],
                                           spr09$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2009',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",line=1.25,side=4)

#2010
spr10<-sprsum[sprsum$Year==2010,] #had most observations
spr10<-spr10[!is.na(spr10$mean_sal),]

zlim<-range(spr10$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr10,span=0.35,degree=2)
lond<-seq(min(spr10$Longitude,na.rm=T),max(spr10$Longitude,na.rm=T),length=100)
latd<-seq(min(spr10$Latitude,na.rm=T),max(spr10$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr10$Latitude[j],
                                           spr10$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2010',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2011
spr11<-sprsum[sprsum$Year==2011,] #had most observations
spr11<-spr11[!is.na(spr11$mean_sal),]

zlim<-range(spr11$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr11,span=0.35,degree=2)
lond<-seq(min(spr11$Longitude,na.rm=T),max(spr11$Longitude,na.rm=T),length=100)
latd<-seq(min(spr11$Latitude,na.rm=T),max(spr11$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr11$Latitude[j],
                                           spr11$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2011',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",line=1.25,side=4)

#2012
spr12<-sprsum[sprsum$Year==2012,] #had most observations
spr12<-spr12[!is.na(spr12$mean_sal),]

zlim<-range(spr12$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr12,span=0.35,degree=2)
lond<-seq(min(spr12$Longitude,na.rm=T),max(spr12$Longitude,na.rm=T),length=100)
latd<-seq(min(spr12$Latitude,na.rm=T),max(spr12$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr12$Latitude[j],
                                           spr12$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2012',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2013
spr13<-sprsum[sprsum$Year==2013,] #had most observations
spr13<-spr13[!is.na(spr13$mean_sal),]

zlim<-range(spr13$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr13,span=0.35,degree=2)
lond<-seq(min(spr13$Longitude,na.rm=T),max(spr13$Longitude,na.rm=T),length=100)
latd<-seq(min(spr13$Latitude,na.rm=T),max(spr13$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr13$Latitude[j],
                                           spr13$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2013',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)


#2014
spr14<-sprsum[sprsum$Year==2014,] 
spr14<-spr14[!is.na(spr14$mean_sal),]

zlim<-range(spr14$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr14,span=0.35,degree=2)
lond<-seq(min(spr14$Longitude,na.rm=T),max(spr14$Longitude,na.rm=T),length=100)
latd<-seq(min(spr14$Latitude,na.rm=T),max(spr14$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr14$Latitude[j],
                                           spr14$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=18,height=15)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2014',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2015
spr15<-sprsum[sprsum$Year==2015,] #had most observations
spr15<-spr15[!is.na(spr15$mean_sal),]

zlim<-range(spr15$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr15,span=0.35,degree=2)
lond<-seq(min(spr15$Longitude,na.rm=T),max(spr15$Longitude,na.rm=T),length=100)
latd<-seq(min(spr15$Latitude,na.rm=T),max(spr15$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr15$Latitude[j],
                                           spr15$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2015',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2016
spr16<-sprsum[sprsum$Year==2016,] 
spr16<-spr16[!is.na(spr16$mean_sal),]

zlim<-range(spr16$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr16,span=0.35,degree=2)
lond<-seq(min(spr16$Longitude,na.rm=T),max(spr16$Longitude,na.rm=T),length=100)
latd<-seq(min(spr16$Latitude,na.rm=T),max(spr16$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr16$Latitude[j],
                                           spr16$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2016',
      xlim=c(-180,-155),ylim=c(52,63))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",side=4,line=1.25)

#2017
spr17<-sprsum[sprsum$Year==2017,] #had most observations
spr17<-spr17[!is.na(spr17$mean_sal),]

zlim<-range(spr17$mean_sal,na.rm=T)

loess<-loess(mean_sal~Longitude*Latitude,data=spr17,span=0.35,degree=2)
lond<-seq(min(spr17$Longitude,na.rm=T),max(spr17$Longitude,na.rm=T),length=100)
latd<-seq(min(spr17$Latitude,na.rm=T),max(spr17$Latitude,na.rm=T),length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist<-min(distance.function(predict.grid$Latitude[j],
                                           predict.grid$Longitude[j],
                                           spr17$Latitude[j],
                                           spr17$Longitude[j]))}

sss.pred<-predict(loess,newdata=predict.grid)

windows(width=6,height=5)
par(oma=c(1,1,1,4))
image(lond,latd,sss.pred,col=hcl.colors(100,"SunsetDark",rev=T),zlim=zlim,
      ylab=expression(paste("Latitude ("^0,'N)')),
      xlab=expression(paste("Longitude ("^0,'E)')),main='2017',
      xlim=c(-180,-155),ylim=c(60,70))
contour(lond,latd,sss.pred,levels=c(32),add=T,lwd=2)
map("worldHires",fill=T,col="gainsboro",add=T)
par(oma=c(1,1,1,1.5))
image.plot(legend.only=T,zlim=zlim,col=hcl.colors(100,"SunsetDark",rev=T),
           legend.shrink=0.3,legend.width=1)
mtext("SSS (psu)",line=1.25,side=4)











