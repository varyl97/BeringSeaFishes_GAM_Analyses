## Script for some figures and tables: 

# Simple figure for domain annotations 
plot(1,1,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),xlim=c(-180,-155),ylim=c(52,65))
abline(h=60,col="firebrick4",lwd=2)
contour(bathy,levels=-c(50,100),labcex=0.8,col='black',add=T)
map("world",fill=T,col="gainsboro",add=T)
raster::scalebar(d=100,xy=c(-158,52),type="bar",divs=2,lonlat=TRUE,label=c(0,NA,100),lwd=3,adj=c(0,-0.75),cex=0.6,below="km")

# for plotting current trajectories over it
plot(1,1,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),xlim=c(-180,-155),ylim=c(52,65))
contour(bathy,levels=-c(50,100),labcex=0.8,col='grey45',add=T)
map("world",fill=T,col="gainsboro",add=T)
raster::scalebar(d=100,xy=c(-158,52),type="bar",divs=2,lonlat=TRUE,label=c(0,NA,100),lwd=3,adj=c(0,-0.75),cex=0.6,below="km")

# Visualization of sampling: 
#for eggs first

pksub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkEggs.csv',header=TRUE,
                check.names=TRUE)
pksub<-pksub[!is.na(pksub$lon)&!is.na(pksub$lat),]

nlat=30;nlon=25
latd=seq(min(pksub$lat,na.rm=T),max(pksub$lat,na.rm=T),length.out=nlat)
lond=seq(min(pksub$lon,na.rm=T),max(pksub$lon,na.rm=T),length.out=nlat)

grid.lon=data.frame(
  lon1=rep(lond[-length(lond)],(nlat-1)),
  lon2=rep(lond[-1],(nlat-1)),
  lon3=rep(lond[-1],(nlat-1)),
  lon4=rep(lond[-length(lond)],(nlat-1)))

grid.lat=data.frame(
  lat1=sort(rep(latd[-length(latd)],(nlon-1))),
  lat2=sort(rep(latd[-length(latd)],(nlon-1))),
  lat3=sort(rep(latd[-1],(nlon-1))),
  lat4=sort(rep(latd[-1],(nlon-1))))

n.stations=NA*(1:nrow(grid.lon))
n.years=NA*(1:nrow(grid.lon))

windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
plot(1,1,ylim=c(50,62),xlim=c(-180,-155),ylab=expression(paste("Latitude ("^0,'N)')),
     xlab=expression(paste("Longitude ("^0,'E)')))

for(i in 1:length(n.stations)){
  tmp=in.chull(pksub$lon,pksub$lat,grid.lon[i,],grid.lat[i,])
  n.years[i]=length(unique(pksub$year[tmp]))
  n.stations[i]=length(pksub$year[tmp])
  points(pksub$lon[tmp],pksub$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #not done yet

z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
z.matrix.stations<-matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F)  
z.matrix.years<-matrix(n.years,ncol=length(z.lat),nrow=length(z.lon),byrow=F)


windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
image.plot(z.lon,z.lat,z.matrix.stations,col=viridis(30),
           xlim=c(-180,-155),ylim=c(52,63), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) # number of stations sampled

windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
image.plot(z.lon,z.lat,z.matrix.years,col=viridis(22),
           xlim=c(-180,-155),ylim=c(52,63), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) # number of years sampled

#for larvae
pklarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',
                     header=TRUE,check.names=TRUE)

nlat=30;nlon=25
latd=seq(min(pklarv.ctd$lat,na.rm=T),max(pklarv.ctd$lat,na.rm=T),length.out=nlat)
lond=seq(min(pklarv.ctd$lon,na.rm=T),max(pklarv.ctd$lon,na.rm=T),length.out=nlat)

grid.lon=data.frame(
  lon1=rep(lond[-length(lond)],(nlat-1)),
  lon2=rep(lond[-1],(nlat-1)),
  lon3=rep(lond[-1],(nlat-1)),
  lon4=rep(lond[-length(lond)],(nlat-1)))

grid.lat=data.frame(
  lat1=sort(rep(latd[-length(latd)],(nlon-1))),
  lat2=sort(rep(latd[-length(latd)],(nlon-1))),
  lat3=sort(rep(latd[-1],(nlon-1))),
  lat4=sort(rep(latd[-1],(nlon-1))))

n.stations=NA*(1:nrow(grid.lon))
n.years=NA*(1:nrow(grid.lon))

windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
plot(1,1,ylim=c(50,62),xlim=c(-180,-155),ylab=expression(paste("Latitude ("^0,'N)')),
     xlab=expression(paste("Longitude ("^0,'E)')))

for(i in 1:length(n.stations)){
  tmp=in.chull(pklarv.ctd$lon,pklarv.ctd$lat,grid.lon[i,],grid.lat[i,])
  n.years[i]=length(unique(pklarv.ctd$year[tmp]))
  n.stations[i]=length(pklarv.ctd$year[tmp])
  points(pklarv.ctd$lon[tmp],pklarv.ctd$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #not done yet

z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
z.matrix.stations<-matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F)  
z.matrix.years<-matrix(n.years,ncol=length(z.lat),nrow=length(z.lon),byrow=F)


windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
image.plot(z.lon,z.lat,z.matrix.stations,col=viridis(30),
           xlim=c(-180,-155),ylim=c(52,63), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #number of stations sampled

windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
image.plot(z.lon,z.lat,z.matrix.years,col=viridis(22),
           xlim=c(-180,-155),ylim=c(52,63), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) # number of years sampled

# Spatial coverage of reg.sst 
mar.sst<-read.csv('./Environmental Data/Mar_SST_ByLocation_NCEP_BS.csv')

windows(width=6,height=5)
par(mai=c(1,1,0.5,0.9))
plot(1,1,xlim=c(-180,-155),ylim=c(51,62),
     ylab=expression(paste("Latitude ("^0,'N)')),
     xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
points(mar.sst$lon,mar.sst$lat,col='grey15',pch=16)
map("worldHires",fill=T,col="gainsboro",add=T)

#Cruise information 

apsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApEggs.csv',header=TRUE,
                              check.names=TRUE)

tmp<-sort(unique(apsub$CRUISE))

for(i in 1:length(tmp)){
  print(tmp[i])
  print(range(apsub$GMT_DATE_TIME[apsub$CRUISE==tmp[i]]))
  print(dim(apsub[apsub$CRUISE==tmp[i],]))
}

# Boxplots for species' breadths: 
ak_clip<-read.csv('./Ichthyo Data/Ak_larv_clip_60.csv',header=TRUE,check.names=TRUE)
ak_clip$species<-'AP'
pk_clip<-read.csv('./Ichthyo Data/Pk_larv_clip_60.csv',header=TRUE,check.names=TRUE)
pk_clip$species<-'WP'
pc_clip<-read.csv('./Ichthyo Data/Pc_larv_clip_60.csv',header=TRUE,check.names=TRUE)
pc_clip$species<-'PC'
yf_clip<-read.csv('./Ichthyo Data/Yfs_larv_clip_60.csv',header=TRUE,check.names=TRUE)
yf_clip$species<-'YFS'
nr_clip<-read.csv('./Ichthyo Data/Nrs_larv_clip_60.csv',header=TRUE,check.names=TRUE)
nr_clip$species<-'NRS'
fh_clip<-read.csv('./Ichthyo Data/Fh_larv_clip_60.csv',header=TRUE,check.names=TRUE)
fh_clip$species<-'FHS'

all_clip<-rbind(ak_clip,pk_clip,pc_clip,yf_clip,nr_clip,fh_clip)

temp_all<-ggplot(all_clip,aes(x=species,y=temperature,group=species))+geom_boxplot()+   
  theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1.2))+
  labs(title="Predictions against Temperature",x='Species',y='Temperature')+
  scale_y_continuous(breaks=round(seq(min(all_clip$temperature),max(all_clip$temperature),by=0.75),1))

sal_all<-ggplot(all_clip,aes(x=species,y=salinity,group=species))+geom_boxplot()+
  theme_bw()+
  labs(title="Predictions against Salinity",x='Species',y='Salinity')+
  scale_y_continuous(breaks=round(seq(min(all_clip$salinity),max(all_clip$salinity),by=0.5),1))

windows(width=10,height=8)
par(mai=c(1,1,0.5,0.9))
ggarrange(temp_all,sal_all)






