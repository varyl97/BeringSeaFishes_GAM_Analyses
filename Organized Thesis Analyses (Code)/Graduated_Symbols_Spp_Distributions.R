## The following code is to create a map of egg/larvae distributions with 
  # graduated symbols that have a legend which matches the graduated 
  # levels automatically (rather than manually, as have done previously). 

## I prefer base R graphics, but seems ggplot is best way to go for this: 

# BS map for ggplot()
world<-map_data("world2Hires")
BSmap<-world[world$long>180&world$long<210&world$lat>50&world$lat<63,]
BSmap$long<-BSmap$long-360
#geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))

appos<-apsub[apsub$Cper10m2>0,]
appos$logcpue<-log(appos$Cper10m2)
apzer<-apsub[apsub$Cper10m2==0,]

windows()
#ggplot(BSmap,aes(x=long,y=lat))+
 # coord_fixed(xlim=c(-180,-150),ylim=c(50,63),ratio=1.3)+
  #geom_polygon(fill="gainsboro",color="darkgrey")+
  #geom_point(data=appos,aes(x=lon,y=lat,size=logcpue),
  #                  fill="steelblue4",shape=21,alpha=0.8)+
  #scale_size_continuous(range=c(0,9),breaks=pretty_breaks(9))+
  #theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
   #     panel.background = element_blank(), axis.line = element_line(colour = "black"))


library(ggmap)
albox<-make_bbox(lon=apsub$lon,lat=apsub$lat,f=0.1)

almap<-get_stamenmap(bbox=albox,maptype="toner-lite",source="google") #this is cool
almap<-readRDS("./Environmental Data/Alaska_Grey_Map.rds") #save a LOT of time

windows()
ggmap(almap)+geom_point(data=appos,aes(x=lon,y=lat,size=logcpue),
                        fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous(range=c(0,9),breaks=pretty_breaks(9))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


apysub<-apsub%>%mutate(yearbin=case_when(year<=1984~'1979-1984',
                                         year<=1989&year>1984~'1985-1989',
                                         year<=1994&year>1989~'1990-1994',
                                         year<=1999&year>1994~'1995-1999',
                                         year<=2004&year>1999~'2000-2004',
                                         year<=2009&year>2004~'2005-2009',
                                         year<=2014&year>2009~'2010-2014',
                                         year<=2019&year>2014~'2015-2016'))
apypos<-apysub[apysub$Cper10m2>0,]
apyz<-apysub[apysub$Cper10m2==0,]
apypos$logcpue<-log(apypos$Cper10m2)

tmp1<-sort(unique(apysub$yearbin))

apy1<-apypos[apypos$yearbin==tmp1[[1]],]
apz1<-apyz[apyz$yearbin==tmp1[[1]],]

apy2<-apypos[apypos$yearbin==tmp1[[2]],]
apz2<-apyz[apyz$yearbin==tmp1[[2]],]

apy3<-apypos[apypos$yearbin==tmp1[[3]],]
apz3<-apyz[apyz$yearbin==tmp1[[3]],]

apy4<-apypos[apypos$yearbin==tmp1[[4]],]
apz4<-apyz[apyz$yearbin==tmp1[[4]],]

apy5<-apypos[apypos$yearbin==tmp1[[5]],]
apz5<-apyz[apyz$yearbin==tmp1[[5]],]

apy6<-apypos[apypos$yearbin==tmp1[[6]],]
apz6<-apyz[apyz$yearbin==tmp1[[6]],]

apy7<-apypos[apypos$yearbin==tmp1[[7]],]
apz7<-apyz[apyz$yearbin==tmp1[[7]],]

apy8<-apypos[apypos$yearbin==tmp1[[8]],]
apz8<-apyz[apyz$yearbin==tmp1[[8]],]


bathybath<-as.bathy(bathy)
  
windows()
year1<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth (m)",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy1,aes(x=lon,y=lat,size=logcpue),
                        fill="#993300",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz1,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  labs(title="Plaice Eggs, 1979-1984",shape="Log CPUE")+ylab("Latitude")+xlab("Longitude")

year2<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy2,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz2,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 1985-1989")+ylab("Latitude")+xlab("Longitude")

year3<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy3,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz3,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 1990-1994")+ylab("Latitude")+xlab("Longitude")

year4<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy4,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz4,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 1995-1999")+ylab("Latitude")+xlab("Longitude")

year5<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy5,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz5,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 2000-2004")+ylab("Latitude")+xlab("Longitude")

year6<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy6,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz6,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 2005-2009")+ylab("Latitude")+xlab("Longitude")

year7<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy7,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz7,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 2010-2014")+ylab("Latitude")+xlab("Longitude")

year8<-autoplot.bathy(bathybath,geom=c("contour","raster"),coast=TRUE,color="grey")+
  scale_fill_gradientn("Depth",values=scales::rescale(c(-7600,0,1,2321)),
                       colors=c("steelblue4","#C7E0FF","grey50","grey80"))+
  geom_point(data=apy8,aes(x=lon,y=lat,size=logcpue),
                               fill="steelblue4",shape=21,alpha=0.8)+
  scale_size_continuous("Log CPUE",range=c(0,9),breaks=pretty_breaks(9))+
  geom_point(data=apz8,aes(x=lon,y=lat),fill="violetred4",shape=8)+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=8))+
  scale_shape_manual(values=('No Catch'=8))+
  labs(title="Plaice Eggs, 2015-2016")+ylab("Latitude")+xlab("Longitude")

windows(width=15,height=15)
year1

windows(width=15,height=15)
ggarrange(year2,year3,year4,nrow=2,ncol=2,common.legend=TRUE,legend="bottom")

windows(width=15,height=15)
ggarrange(year5,year6,year7,year8,ncol=2,nrow=2,common.legend=TRUE,legend="bottom")

years<-sort(unique(apypos$yearbin))
tmp1<-1:ceiling(length(years)/4)
for(i in 1:length(apypos)){
  par(omi=c(0.25,0.3,0.55,0.25))
      apeplots=ggmap(almap)+geom_point(data=apypos%>%filter(yearbin==years[i]),
                                  aes(x=lon,y=lat,
                                      size=logcpue%>%filter(yearbin==years[i])),
                            fill="steelblue4",shape=21,alpha=0.8)+
      theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))+
      labs(title=as.character(years[i]))+
      xlab("Longitude")+ylab("Latitude")
  scale_size_continuous(range=c(0,8),breaks=pretty_breaks(8))
}


for(j in 1:length(tmp1)){
  par(mfrow=c(2,2),omi=c(0.25,0.3,0.55,0.25),mar=c(2,2,1.5,1.5))
  for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
    plot(1,1,xlim=c(-180,-155),ylim=c(52,63),
         main=as.character(years[i]),ylab=expression(paste("Latitude ("^0,'N)')),
         xlab=expression(paste("Longitude ("^0,'E)')))
    contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
    map("worldHires",fill=T,col="gainsboro",add=T)
    symbols(apysub$lon[apysub$Cper10m2>0&apysub$yearbin==years[i]],
            apysub$lat[apysub$Cper10m2>0&apysub$yearbin==years[i]],
            circles=apysub$Cper10m2[apysub$Cper10m2>0&apysub$yearbin==years[i]],
            inches=0.25*max(apysub$Cper10m2[apysub$Cper10m2>0&apysub$yearbin==years[i]])/max_density,
            bg="steelblue2",add=T)
    points(apysub$lon[apysub$Cper10m2==0&apysub$yearbin==years[i]],
           apysub$lat[apysub$Cper10m2==0&apysub$yearbin==years[i]],
           pch="+",col="violetred4",cex=0.4)
    mtext("Plaice Egg Catches", outer=T,cex=1,line=1)
  }
  legend("bottomright",legend=lg.text,pch=c(3,16,16,16,16),pt.cex=c(0.8,0.8,1.4,2,4),
         col=c("violetred4",rep("steelblue2",4)))}
