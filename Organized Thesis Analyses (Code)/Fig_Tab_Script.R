## Script for some figures and tables: 

# Simple figure for domain annotations 
plot(1,1,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),xlim=c(-180,-155),ylim=c(52,65))
abline(h=60,col="firebrick4",lwd=2)
contour(bathy,levels=-c(50,100),labcex=0.8,col='black',add=T)
map("world",fill=T,col="gainsboro",add=T)
raster::scalebar(d=100,xy=c(-158,52),type="bar",divs=2,lonlat=TRUE,label=c(0,NA,100),lwd=3,adj=c(0,-0.75),cex=0.6,below="km")

# for plotting current trajectories over it
plot(1,1,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),xlim=c(-180,-155),ylim=c(52,65))
contour(bathy,levels=-c(50,100,200),labcex=0.75,col='grey45',add=T)
map("world",fill=T,col="gainsboro",add=T)
raster::scalebar(d=100,xy=c(-158,52),type="bar",divs=2,lonlat=TRUE,label=c(0,NA,100),lwd=3,adj=c(0,-0.75),cex=0.6,below="km")

# for plotting regional temperature index
par(mai=c(1,1,0.5,0.9))
plot(1,1,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),xlim=c(-180,-155),ylim=c(52,65))
contour(bathy,levels=-c(50,100,200),labcex=0.75,col='grey70',add=T)
points(mar.sst.loc$lon,mar.sst.loc$lat,pch=18,cex=1.5,col="dodgerblue3")
map("world",fill=T,col="gainsboro",add=T)
raster::scalebar(d=100,xy=c(-158,52),type="bar",divs=2,lonlat=TRUE,label=c(0,NA,100),lwd=3,adj=c(0,-0.75),cex=0.6,below="km")

# Visualization of sampling: 
#for eggs first

pksub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkEggs.csv',header=TRUE,
                check.names=TRUE)
pksub<-pksub[!is.na(pksub$lon)&!is.na(pksub$lat),]

nlat=30;nlon=25
latd=seq(min(pkegg$LAT,na.rm=T),max(pkegg$LAT,na.rm=T),length.out=nlat)
lond=seq(min(pkegg$LON,na.rm=T),max(pkegg$LON,na.rm=T),length.out=nlat)

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
  tmp=in.chull(pkegg$LON,pkegg$LAT,grid.lon[i,],grid.lat[i,])
  n.years[i]=length(unique(pkegg$YEAR_[tmp]))
  n.stations[i]=length(pkegg$YEAR_[tmp])
  points(pkegg$LON[tmp],pkegg$LAT[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #not done yet

z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
z.matrix.stations<-matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F)  
z.matrix.years<-matrix(n.years,ncol=length(z.lat),nrow=length(z.lon),byrow=F)


windows(width=10,height=4);par(mai=c(1,1,0.5,1.5),mfrow=c(1,2))
image.plot(z.lon,z.lat,z.matrix.stations,col=viridis(30),
           xlim=c(-174.5,-155),ylim=c(53,60), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),
           legend.lab='S')
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) # number of stations sampled

image.plot(z.lon,z.lat,z.matrix.years,col=viridis(22),
           xlim=c(-174.5,-155),ylim=c(53,60), main="",
           ylab="",
           xlab=expression(paste("Longitude ("^0,'E)')),
           legend.lab="Y")
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) # number of years sampled

#for larvae
pklarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',
                     header=TRUE,check.names=TRUE)

nlat=30;nlon=25
latd=seq(min(yflarv.ctd$lat,na.rm=T),max(yflarv.ctd$lat,na.rm=T),length.out=nlat)
lond=seq(min(yflarv.ctd$lon,na.rm=T),max(yflarv.ctd$lon,na.rm=T),length.out=nlat)

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
  tmp=in.chull(yflarv.ctd$lon,yflarv.ctd$lat,grid.lon[i,],grid.lat[i,])
  n.years[i]=length(unique(yflarv.ctd$year[tmp]))
  n.stations[i]=length(yflarv.ctd$year[tmp])
  points(yflarv.ctd$lon[tmp],yflarv.ctd$lat[tmp],col=i,pch=16)
  polygon(grid.lon[i,],grid.lat[i,])
}
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #not done yet

z.lat<-(latd[1:(length(latd)-1)]+latd[2:length(latd)])/2
z.lon<-(lond[1:(length(lond)-1)]+lond[2:length(lond)])/2
z.matrix.stations<-matrix(n.stations,ncol=length(z.lat),nrow=length(z.lon),byrow=F)  
z.matrix.years<-matrix(n.years,ncol=length(z.lat),nrow=length(z.lon),byrow=F)


windows(width=15,height=5);par(mai=c(1,1,0.5,1.5),mfrow=c(1,2))
image.plot(z.lon,z.lat,z.matrix.stations,col=viridis(30),
           xlim=c(-180,-155),ylim=c(52,65), main="",
           ylab=expression(paste("Latitude ("^0,'N)')),
           legend.lab="Stations Sampled",
           xlab=expression(paste("Longitude ("^0,'E)')))
contour(bathy,levels=-c(100,200),cex=0.5,col='grey28',add=T)
map("worldHires",fill=T,col="gainsboro",add=T) #number of stations sampled

image.plot(z.lon,z.lat,z.matrix.years,col=viridis(22),
           xlim=c(-180,-155),ylim=c(52,65), main="",
           ylab="",legend.lab="Years Sampled",
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

#Cruises 2.0 - simplified: 
tmp<-sort(unique(apsub$year))
for(i in 1:length(tmp)){
  print(sort(unique(apsub$CRUISE[apsub$year==tmp[i]])))
  print(range(apsub$GMT_DATE_TIME[apsub$year==tmp[i]]))
  print(unique(apsub$SS[apsub$year==tmp[i]]))
}

tmp<-sort(unique(aplarv.ctd$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(aplarv.ctd$CRUISE[aplarv.ctd$year==tmp[i]])))
  print(range(aplarv.ctd$GMT_DATE_TIME[aplarv.ctd$year==tmp[i]]))
  print(unique(aplarv.ctd$SS[aplarv.ctd$year==tmp[i]]))
}

tmp<-sort(unique(yflarv.ctd$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(yflarv.ctd$CRUISE[yflarv.ctd$year==tmp[i]])))
  print(range(yflarv.ctd$GMT_DATE_TIME[yflarv.ctd$year==tmp[i]]))
  print(unique(yflarv.ctd$SS[yflarv.ctd$year==tmp[i]]))
}

tmp<-sort(unique(fhsub$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(fhsub$CRUISE[fhsub$year==tmp[i]])))
  print(range(fhsub$GMT_DATE_TIME[fhsub$year==tmp[i]]))
  print(unique(fhsub$SS[fhsub$year==tmp[i]]))
}

tmp<-sort(unique(fhlarv.ctd$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(fhlarv.ctd$CRUISE[fhlarv.ctd$year==tmp[i]])))
  print(range(fhlarv.ctd$date[fhlarv.ctd$year==tmp[i]]))
  print(unique(fhlarv.ctd$SS[fhlarv.ctd$year==tmp[i]]))
}

tmp<-sort(unique(pksub$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(pksub$CRUISE[pksub$year==tmp[i]])))
  print(range(pksub$date[pksub$year==tmp[i]]))
  print(unique(pksub$SS[pksub$year==tmp[i]]))
}

tmp<-sort(unique(pklarv.ctd$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(pklarv.ctd$CRUISE[pklarv.ctd$year==tmp[i]])))
  print(range(pklarv.ctd$date[pklarv.ctd$year==tmp[i]]))
  print(unique(pklarv.ctd$SS[pklarv.ctd$year==tmp[i]]))
}

tmp<-sort(unique(rxsub$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(rxsub$CRUISE[rxsub$year==tmp[i]])))
  print(range(rxsub$date[rxsub$year==tmp[i]]))
  print(unique(rxsub$SS[rxsub$year==tmp[i]]))
}

#both nrs & pc 
nrspc<-rbind(nrslarv.ctd,pclarv.ctd)

tmp<-sort(unique(nrspc$year))
for(i in 1:length(tmp)){
  print(tmp[i])
  print(sort(unique(nrspc$CRUISE[nrspc$year==tmp[i]])))
  print(range(nrspc$date[nrspc$year==tmp[i]]))
  print(unique(nrspc$SS[nrspc$year==tmp[i]]))
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


##Scatter plot of temperature and salinity by depth: 
#untrimmed fhlarv.ctd

#warm regimes: 2001 thru 2005 and 2014 thru 2016
#cold regimes: 2006 thru 2013 (from Baker 2021 & Stabeno et al. 2017) 
ctdsub<-pklarv.ctd[pklarv.ctd$month==5,]

ctdsub<-ctdsub %>% mutate(domain = 
                    case_when(bottom_depth>100 ~ "Outer",
                              bottom_depth<=100 & bottom_depth>50 ~ "Middle",
                              bottom_depth<=50 ~ "Inner"))
table(ctdsub$domain) #check 

windows();par(mai=c(1,1,0.5,0.9))
plot(1,1,xlim=c(29,33),ylim=c(-1.7,14),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Inner"],
       ctdsub$temperature[ctdsub$domain=="Inner"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Middle"],
       ctdsub$temperature[ctdsub$domain=="Middle"],
       pch=16,col="#2c92e6")
points(ctdsub$salinity[ctdsub$domain=="Outer"],
       ctdsub$temperature[ctdsub$domain=="Outer"],
       pch=16,col="#054375")

ctdsub<-ctdsub %>% mutate(regime = 
                            case_when(year>=2006 & year<=2013 ~ "Cold",
                                      year>=2001 & year<=2005 ~ "Warm1",
                                      year>=2014 & year<=2016 ~ "Warm2"))
table(ctdsub$regime)

windows(width=10,height=5);par(oma=c(0.5,2.5,0.5,2.5),mfrow=c(1,3))
plot(1,1,xlim=c(29,33),ylim=c(-1.7,14),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Warm1"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Warm1"],
       pch=16,col="#2c92e6")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Warm1"],
       pch=16,col="#054375")

plot(1,1,xlim=c(29,33),ylim=c(-1.7,14),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Cold"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Cold"],
       pch=16,col="#2c92e6")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Cold"],
       pch=16,col="#054375")

plot(1,1,xlim=c(29,33),ylim=c(-1.7,14),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Warm2"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Warm2"],
       pch=16,col="#2c92e6")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Warm2"],
       pch=16,col="#054375")
legend(x="bottomleft",legend=c("Inner","Middle","Outer"),
       col=c("#9ac2e3","#2c92e6","#054375"),pch=c(16,16,16))

#visualization 2
windows(width=10,height=5);par(oma=c(0.5,1.75,0.5,1.75),mfrow=c(1,3))
plot(1,1,xlim=range(ctdsub$salinity,na.rm=T),
     ylim=range(ctdsub$temperature,na.rm=T),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Warm1"],
       pch=16,col="#750507")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Cold"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Inner"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Inner"&ctdsub$regime=="Warm2"],
       pch=16,col="#750507")

plot(1,1,xlim=range(ctdsub$salinity,na.rm=T),
     ylim=range(ctdsub$temperature,na.rm=T),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Warm1"],
       pch=16,col="#750507")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Cold"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Middle"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Middle"&ctdsub$regime=="Warm2"],
       pch=16,col="#750507")

plot(1,1,xlim=range(ctdsub$salinity,na.rm=T),
     ylim=range(ctdsub$temperature,na.rm=T),
     ylab=expression(paste("Temperature (" ^0,"C)")),
     xlab="Salinity (psu)")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Warm1"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Warm1"],
       pch=16,col="#750507")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Cold"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Cold"],
       pch=16,col="#9ac2e3")
points(ctdsub$salinity[ctdsub$domain=="Outer"&ctdsub$regime=="Warm2"],
       ctdsub$temperature[ctdsub$domain=="Outer"&ctdsub$regime=="Warm2"],
       pch=16,col="#750507")
legend(x="bottomleft",legend=c("Warm Regime","Cold Regime"),
       col=c("#750507","#9ac2e3"),pch=c(16,16))


windows(width=10,height=5);par(oma=c(0.5,2.5,0.5,2.5),mfrow=c(1,3))

## figure of doy breadths for each species & life stage 
apsub$spestg<-as.numeric(1)
apsubpos<-apsub[apsub$Cper10m2>0,]
apsubpos<-apsubpos[c("doy","spestg")]
aplarv.ctd$spestg<-as.numeric(2)
aplarvpos<-aplarv.ctd[aplarv.ctd$Cper10m2>0,]
aplarvpos<-aplarvpos[c("doy","spestg")]

fhsub$spestg<-as.numeric(3)
fhsubpos<-fhsub[fhsub$Cper10m2>0,]
fhsubpos<-fhsubpos[c("doy","spestg")]
fhlarv.ctd$spestg<-as.numeric(4)
fhlarvpos<-fhlarv.ctd[fhlarv.ctd$Cper10m2>0,]
fhlarvpos<-fhlarvpos[c("doy","spestg")]

nrslarv.ctd$spestg<-as.numeric(5)
nrslarvpos<-nrslarv.ctd[nrslarv.ctd$Cper10m2>0,]
nrslarvpos<-nrslarvpos[c("doy","spestg")]

pclarv.ctd$spestg<-as.numeric(6)
pclarvpos<-pclarv.ctd[pclarv.ctd$Cper10m2>0,]
pclarvpos<-pclarvpos[c("doy","spestg")]

pksub$spestg<-as.numeric(7)
pksubpos<-pksub[pksub$Cper10m2>0,]
pksubpos<-pksubpos[c("doy","spestg")]
pklarv.ctd$spestg<-as.numeric(8)
pklarvpos<-pklarv.ctd[pklarv.ctd$Cper10m2>0,]
pklarvpos<-pklarvpos[c("doy","spestg")]

yflarv.ctd$spestg<-as.numeric(9)
yflarvpos<-yflarv.ctd[yflarv.ctd$Cper10m2>0,]
yflarvpos<-yflarvpos[c("doy","spestg")]

rxsub$spestg<-as.numeric(10)
rxsubpos<-rxsub[rxsub$Cper10m2>0,]
rxsubpos<-rxsubpos[c("doy","spestg")]

windows(width=8,height=10);par(mai=c(1,1,0.5,0.9))
plot(1,1,xlab="Day of Year",ylab="",
     xlim=c(90,285),ylim=c(1,10))
points(apsubpos$doy,apsubpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(aplarvpos$doy,aplarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(fhsubpos$doy,fhsubpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(fhlarvpos$doy,fhlarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(nrslarvpos$doy,nrslarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(pclarvpos$doy,pclarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(pksubpos$doy,pksubpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(pklarvpos$doy,pklarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(yflarvpos$doy,yflarvpos$spestg,pch=16,col="#8db6c7",cex=0.8)
points(rxsubpos$doy,rxsubpos$spestg,pch=16,col="#8db6c7",cex=0.8)


