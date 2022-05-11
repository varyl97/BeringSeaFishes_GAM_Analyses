
#Script for all year temp and sal plots 
# 4/28/2022

muctd$month<-str_sub(muctd$Date,start=1,end=1) #not right for double digit months but doesn't matter 
cutctd<-muctd[muctd$month>3&muctd$month<6,]

zlim.t<-range(cutctd$mean_temp,na.rm=T)
zlim.s<-range(cutctd$mean_sal,na.rm=T)

loess.t<-loess(mean_temp~Longitude*Latitude,
               data=cutctd,degree=2,span=0.15)
loess.s<-loess(mean_sal~Longitude*Latitude,
               data=cutctd,degree=2,span=0.15)

lond<-seq(min(cutctd$Longitude,na.rm=T),max(cutctd$Longitude,na.rm=T),
          length=100)
latd<-seq(min(cutctd$Latitude,na.rm=T),max(cutctd$Latitude,na.rm=T),
          length=100)
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c("Longitude","Latitude")

for(j in 1:nrow(predict.grid)){
  predict.grid$dist[j]<-min(distance.function(
    predict.grid$Latitude[j],predict.grid$Longitude[j],
    cutctd$Latitude[j],cutctd$Longitude[j]))}
predict.grid[predict.grid$dist>744438,]<-NA

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
points(cutctd$Longitude,cutctd$Latitude,pch=16)
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
points(cutctd$Longitude,cutctd$Latitude,pch=16)
map("worldHires",fill=T,col="gainsboro",add=T)