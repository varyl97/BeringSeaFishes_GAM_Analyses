# T-S Diagrams with Species Overlaid --------------------------------------
#This code creates a temperature-salinity diagram (T-S Diagram) for the Eastern Bering Sea and overlays 
#   polygons that represent conditions that are associated with each individual species. 
#This code is slightly modified from Hafez Ahmad's guide 
#(https://hafezahmad.medium.com/making-temperature-salinity-diagrams-called-the-t-s-diagram-with-python-and-r-programming-5deec6378a29)
library(marelac)
library(plot3D)
library(DescTools)

#Load in the CTD data from which temperature and salinity values are obtained: 
dat<-read.csv('./Environmental Data/All_CTD_Data_8302021.csv',header=TRUE,check.names=TRUE)
dat<-subset(dat,Depth<11,na.rm=T)
dat<-dat[dat$Temperature<14,]
dat<-dat[dat$Salinity>29&dat$Salinity<36,]

mint=min(dat$Temperature,na.rm=T)
maxt=max(dat$Temperature,na.rm=T)
mins=min(dat$Salinity,na.rm=T)
maxs=max(dat$Salinity,na.rm=T)

temp<-seq(from=mint,to=maxt,length.out=180)
sal<-seq(from=mins,to=maxs,length.out=180)

sigma.c<-outer(sal,temp,FUN=function(S,t)sw_dens(S=S,t=t)-1000)
sigma.c

pkcol<-adjustcolor('#BFA0C6',alpha.f=0.6)


windows()
par(mai=c(1,1,0.5,0.9))
contour2D(x=sal,y=temp,z=sigma.c,lwd=2,main='Temperature-Salinity Diagram',col='black',
          xlab=expression('Salinity (psu)'),ylab=expression('Temperature'*~degree*C*')'))
DrawEllipse(x=31.5,y=6.35,radius.x=0.5,radius.y=1.15,rot=0, border='grey69',col='coral')#x axis = salinity, y axis = temperature
DrawEllipse(x=31.4,y=5,radius.x=0.4,radius.y=0.5,rot=0,border='grey69',col='cyan4')
DrawEllipse(x=32,y=5,radius.x=1,radius.y=1,rot=0,border='grey69',col=pkcol)
DrawEllipse(x=30.75,y=9.75,radius.x=0.75,radius.y=0.75,border='grey69',col='blue4')
DrawEllipse(x=31.75,y=2.5,radius.x=0.75,radius.y=2,border='grey69',col='darkgoldenrod4')
DrawEllipse(x=31,y=6.625,radius.x=0.75,radius.y=0.875,border='grey69',col='darkorchid4')
legend('bottomleft',legend=c('Flathead Sole','Alaska Plaice','Walleye Pollock','Yellowfin Sole',
                              'Northern Rock Sole','Pacific Cod'),col=c('coral','cyan4','#BFA0C6','blue4',
                                                                        'darkgoldenrod4','darkorchid4'),
       bg='snow1',pch=16)
