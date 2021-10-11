setwd("C:/Users/varyl/Desktop/MS Research/RStudio")

library(mgcv)
library(sgeostat)
library(colorRamps) 
library(itsadug)   
library(maps)
library(mapdata)
library(marmap)
library(spacetime)
library(fields)
library(date)
library(chron)
library(marmap)
library(raster)
library(ggplot2)
source('distance.function.R')
source("vis.gam_COLORS.R")
library(RColorBrewer)
library(ggrepel)
library(viridisLite)
library(viridis)
library(scales)
library(ggpubr)
library(gganimate)
library(ggmap)
library(shiny)
library(leaflet)
library(spData)
library(sf)
library(dplyr)
library(tmap)
library(cartogram)
library(devtools)
library(spatialEco)
library(tmaptools)
library(rgdal)
library(RStoolbox)
library(MASS)
library(plotly)
library(reshape2)
library(ncdf4)
library(lattice)
#library(ggOceanMaps)
#library(ggOceanMapsData)
library(colorspace)
library(ggstream)
#library(xlsx)
library(plyr)
library(stringr)
library(lubridate)
library(itsadug)
library(RColorBrewer)

###BATHYMETRY###################################################################
str_name<-"BeringSea_1.tiff"
brs_bathy<-raster(str_name)
brs_bathy2<-as.bathy(brs_bathy)

windows()
plot(brs_bathy2,land=TRUE,deep=c(-5000,-200,0),shallow=c(-1000,-50,0),
     step=c(1000,50,0),lwd=c(0.7,0.7,0.7,0.7),
     lty=c(1,1,1),col=c("grey80","black"),
     image=TRUE,drawlabels=T,
     bpal = list(c(min(brs_bathy2), -500, "purple","blue4"),
                 c(-500, -1, "blue4", "blue", "lightblue"),
                 c(1, max(brs_bathy2), "gray90", "gray10")))

###CTD Loading: ################################################################
allctd<-read.csv(file="All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference

allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,]

###POLLOCK: ####################################################################
pkeggraw<-read.csv(file='BS_PollockEggCatch.csv',header=TRUE,check.names=TRUE)
pklarvraw<-read.csv(file='BS_PollockLarvaeCatch.csv',header=TRUE,check.names=TRUE)

pkeggraw<-pkeggraw[pkeggraw$HAUL_ID!='1SS02 81 1 60BON 2',]
pklarvraw<-pklarvraw[pklarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

pkeggraw$doy<-as.numeric(mdy.date(pkeggraw$MONTH_,pkeggraw$DAY_,1960))
pklarvraw$doy<-as.numeric(mdy.date(pklarvraw$MONTH_,pklarvraw$DAY_,1960))
pkeggraw$id<-paste(pkeggraw$CRUISE,pkeggraw$LAT,pkeggraw$LON,pkeggraw$GMT_DATE_TIME,pkeggraw$MESH,sep="_")
pklarvraw$id<-paste(pklarvraw$CRUISE,pklarvraw$LAT,pklarvraw$LON,pklarvraw$GMT_DATE_TIME,pklarvraw$MESH,sep="_")

pkeggraw[pkeggraw$NET==3,] 
pkegg_1_2<-pkeggraw[pkeggraw$NET<3,]
table(pkegg_1_2$NET)
tmp<-table(pkegg_1_2$id)
id_2<-names(tmp)[tmp==2]
pkegg_1<-pkegg_1_2[pkegg_1_2$NET==1,]
pkegg_2<-pkegg_1_2[pkegg_1_2$id%in%id_2&pkegg_1_2$NET==1,]

pkegg_1net<-pkegg_1_2[!pkegg_1_2$id%in%id_2,]
pkegg<-rbind(pkegg_1net,pkegg_2)
dim(pkegg)#We should end up with n = 6181-27-84= 6070 data points

#net cleaning, larvae: 
pklarvraw[pklarvraw$NET==3,]
pklarv_1_2<-pklarvraw[pklarvraw$NET<3,]
tmp<-table(pklarv_1_2$id)
id_2<-names(tmp)[tmp==2]
pklarv_2<-pklarv_1_2[pklarv_1_2$id%in%id_2&pklarv_1_2$NET==1,]

pklarv_1net<-pklarv_1_2[!pklarv_1_2$id%in%id_2,]
pklarvae<-rbind(pklarv_1net,pklarv_2)
dim(pklarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to Y
pkegg<-pkegg[pkegg$PRIMARY_NET=='Y',]
pklarvae<-pklarvae[pklarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
pkegg$dist<-NA;pklarvae$dist<-NA
for(i in 1:nrow(pkegg)){
  pkegg$dist[i]<-min(distance.function(pkegg$LAT[i],pkegg$LON[i],pkegg$LAT[pkegg$LARVALCATCHPER10M2>0],pkegg$LON[pkegg$LARVALCATCHPER10M2>0]))
  pklarvae$dist[i]<-min(distance.function(pklarvae$LAT[i], pklarvae$LON[i], pklarvae$LAT[pklarvae$LARVALCATCHPER10M2>0], pklarvae$LON[pklarvae$LARVALCATCHPER10M2>0]))
}

attach(pklarvae)
pklarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(pklarvae)

tmp<-table(pklarvae$YEAR_)
index<-match(pklarvae$YEAR_,names(tmp))
pklarvae$SS<-as.numeric(tmp[index])

attach(pkegg)
pkegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(pkegg)

tmp<-table(pkegg$YEAR_)
index<-match(pkegg$YEAR_,names(tmp))
pkegg$SS<-as.numeric(tmp[index])

pkegg$DATE<-NULL
pkegg$DATE<-paste(pkegg$MONTH_,pkegg$DAY_,pkegg$YEAR_,sep="/")

pklarvae$DATE<-NULL
pklarvae$DATE<-paste(pklarvae$MONTH_,pklarvae$DAY_,pklarvae$YEAR_,sep="/")

pksub<-pkegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')]
pksub<-subset(pksub,doy>99&doy<160)
pklarv<-pklarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]

names(pksub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(pklarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

pksub$SSB<-NA
pklarv$SSB<-NA

SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
pkSSB<-SSBdf[SSBdf$Species=='walleye pollock',]
pkSSB<-pkSSB[c('Year','SSB')]

for(i in 1:nrow(pklarv)){
  pklarv$SSB[i]<-pkSSB$SSB[pkSSB$Year==pklarv$year[i]]
}

for(i in 1:nrow(pksub)){
  pksub$SSB[i]<-pkSSB$SSB[pkSSB$Year==pksub$year[i]]
}

pksub.ctd<-pksub
pklarv.ctd<-pklarv

#uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")
#idx.eg<-match(pksub.ctd$CRUISE,uniquecruise$unique_cruise)
#idx.lv<-match(pklarv.ctd$CRUISE,uniquecruise$unique_cruise)
#pksub.ctd$ctd_file_name<-uniquecruise$file_name[idx.eg]
#pklarv.ctd$ctd_file_name<-uniquecruise$file_name[idx.lv]

allctd$Link_ID<-NA
pksub.ctd$Link_ID<-NA
pklarv.ctd$Link_ID<-NA

#allctd$Station<-as.numeric(allctd$Station)#conversion in an attempt to make matching work 
#allctd$Haul<-as.numeric(allctd$Haul)
#pksub.ctd$STATION_NAME<-as.numeric(pksub.ctd$STATION) #warnings introduced for cruises: 2MF92, 10C88, 1MP91, 0MF91, 3MF79, 1DN88, MF862 
#pksub.ctd$HAUL_NAME<-as.numeric(pksub.ctd$HAUL)

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
pksub.ctd$Link_ID<-paste(pksub.ctd$CRUISE,pksub.ctd$STATION,pksub.ctd$HAUL,sep="_")
pklarv.ctd$Link_ID<-paste(pklarv.ctd$CRUISE,pklarv.ctd$STATION,pklarv.ctd$HAUL,sep="_")

#pksub.ctd<-pksub.ctd[!is.na(pksub.ctd$ctd_file_name),] 
pksub.ctd$temperature<-NA
pksub.ctd$salinity<-NA
pksub.ctd$CTD_date<-NA
pksub.ctd$CTD_link<-NA

for(i in 1:nrow(pksub.ctd)){
  tryCatch({
    idx1.eg<-allctd$Link_ID==pksub.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.eg,]
    pksub.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    pksub.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    pksub.ctd$CTD_date[i]<-as.character(ctdcast$Date)
    pksub.ctd$CTD_link[i]<-as.character(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(pksub.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
pksub.ctd<-pksub.ctd[!pksub.ctd$temperature=="NaN",] #dim: 1640 24

#pklarv.ctd<-pklarv.ctd[!is.na(pklarv.ctd$ctd_file_name),] 
pklarv.ctd$temperature<-NA
pklarv.ctd$salinity<-NA
pklarv.ctd$CTD_date<-NA
pklarv.ctd$CTD_link<-NA
pklarv.ctd$CTD_time<-NA

for(i in 1:nrow(pklarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==pklarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    pklarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    pklarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    pklarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    pklarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    pklarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(pklarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
pklarv.ctd<-pklarv.ctd[!pklarv.ctd$temperature=="NaN",]#dim: 2646 24

pklarv.ctd$date_diff<-NA
pklarv.ctd$date<-parse_date_time(pklarv.ctd$date,orders="mdy")
pklarv.ctd$CTD_date<-parse_date_time(pklarv.ctd$CTD_date,orders="mdy")
pklarv.ctd$date_diff<-difftime(pklarv.ctd$date,pklarv.ctd$CTD_date,units="hour")
pklarv.ctd<-pklarv.ctd[which(pklarv.ctd$date_diff>(-6)&pklarv.ctd$date_diff<6),]
dim(pklarv.ctd) #2565

###to avoid doing all the above gymnastics: 
pksub<-read.csv(file='Cleaned_Cut_PkEggs.csv',header=TRUE,check.names=TRUE)
pksub.ctd<-read.csv(file='Cleaned_Cut_PkEggs_wCTD.csv',header=TRUE,check.names=TRUE)
pklarv<-read.csv(file='Cleaned_Cut_PkLarv.csv',header=TRUE,check.names=TRUE)
pklarv.ctd<-read.csv(file='Cleaned_Cut_PkLarv_wCTD.csv',header=TRUE,check.names=TRUE)

#YELLOWFIN: 
yfeggraw<-read.csv(file='YFSole_Egg_Catch.csv',header=TRUE,check.names=TRUE)
yflarvraw<-read.csv(file='YFSole_Larvae_Catch.csv',header=TRUE,check.names=TRUE)

yfeggraw<-yfeggraw[yfeggraw$HAUL_ID!='1SS02 81 1 60BON 2',]
yflarvraw<-yflarvraw[yflarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

yfeggraw$doy<-as.numeric(mdy.date(yfeggraw$MONTH_,yfeggraw$DAY_,1960))
yflarvraw$doy<-as.numeric(mdy.date(yflarvraw$MONTH_,yflarvraw$DAY_,1960))
yfeggraw$id<-paste(yfeggraw$CRUISE,yfeggraw$LAT,yfeggraw$LON,yfeggraw$GMT_DATE_TIME,yfeggraw$MESH,sep="_")
yflarvraw$id<-paste(yflarvraw$CRUISE,yflarvraw$LAT,yflarvraw$LON,yflarvraw$GMT_DATE_TIME,yflarvraw$MESH,sep="_")

yfeggraw[yfeggraw$NET==3,]
yfegg_1_2<-yfeggraw[yfeggraw$NET<3,]
tmp<-table(yfegg_1_2$id)
id_2<-names(tmp)[tmp==2]
yfegg_1<-yfegg_1_2[yfegg_1_2$NET==1,]
yfegg_2<-yfegg_1_2[yfegg_1_2$id%in%id_2&yfegg_1_2$NET==1,]
yfegg_1net<-yfegg_1_2[!yfegg_1_2$id%in%id_2,]
yfegg<-rbind(yfegg_1net,yfegg_2)
dim(yfegg)#We should end up with n = 6181-27-84= 6070 data points

yflarvraw[yflarvraw$NET==3,]
yflarv_1_2<-yflarvraw[yflarvraw$NET<3,]
tmp<-table(yflarv_1_2$id)
id_2<-names(tmp)[tmp==2]
yflarv_2<-yflarv_1_2[yflarv_1_2$id%in%id_2&yflarv_1_2$NET==1,]
yflarv_1net<-yflarv_1_2[!yflarv_1_2$id%in%id_2,]
yflarvae<-rbind(yflarv_1net,yflarv_2)
dim(yflarvae)#We should end up with n = 6181-27-84= 6070 data points

yfegg<-yfegg[yfegg$PRIMARY_NET=='Y',]
yflarvae<-yflarvae[yflarvae$PRIMARY_NET=='Y',]

yfegg$dist<-NA;yflarvae$dist<-NA
for(i in 1:nrow(yfegg)){
  yfegg$dist[i]<-min(distance.function(yfegg$LAT[i],yfegg$LON[i],yfegg$LAT[yfegg$LARVALCATCHPER10M2>0],yfegg$LON[yfegg$LARVALCATCHPER10M2>0]))
  yflarvae$dist[i]<-min(distance.function(yflarvae$LAT[i], yflarvae$LON[i], yflarvae$LAT[yflarvae$LARVALCATCHPER10M2>0], yflarvae$LON[yflarvae$LARVALCATCHPER10M2>0]))
}

attach(yflarvae)
yflarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(yflarvae)

tmp<-table(yflarvae$YEAR_)
index<-match(yflarvae$YEAR_,names(tmp))
yflarvae$SS<-as.numeric(tmp[index])

attach(yfegg)
yfegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(yfegg)

tmp<-table(yfegg$YEAR_)
index<-match(yfegg$YEAR_,names(tmp))
yfegg$SS<-as.numeric(tmp[index])

yfsub<-yfegg[yfegg$MONTH_>3&yfegg$MONTH_<10,]
yfsub<-yfsub[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]
yflarv<-yflarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]

names(yfsub)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')
names(yflarv)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                 'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')

yfsub$SSB<-NA
yflarv$SSB<-NA

SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
yfSSB<-SSBdf[SSBdf$Species=='yellowfin sole',]
yfSSB<-yfSSB[c('Year','SSB')]

for(i in 1:nrow(yflarv)){
  yflarv$SSB[i]<-yfSSB$SSB[yfSSB$Year==yflarv$year[i]]
}

for(i in 1:nrow(yfsub)){
  yfsub$SSB[i]<-yfSSB$SSB[yfSSB$Year==yfsub$year[i]]
}

uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")

idx2<-match(yfsub$CRUISE,uniquecruise$unique_cruise)
yfsub$ctd_file_name<-uniquecruise$file_name[idx2]

yfsub<-na.exclude(yfsub)
yfsub$temperature<-NA  
yfsub$salinity<-NA

for (i in 1:nrow(yfsub)){
  tryCatch({
    filename<-paste0(yfsub$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==yfsub$HAUL_NAME[i]&ctd.file$Station==yfsub$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    yfsub$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    yfsub$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(yfsub$CRUISE[i]),conditionMessage(e), "\n")})
}

idx3<-match(yflarv$CRUISE,uniquecruise$unique_cruise)
yflarv$ctd_file_name<-uniquecruise$file_name[idx3]

yflarv<-na.exclude(yflarv)
yflarv$temperature<-NA  
yflarv$salinity<-NA

for (i in 1:nrow(yflarv)){
  tryCatch({
    filename<-paste0(yflarv$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==yflarv$HAUL_NAME[i]&ctd.file$Station==yflarv$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    yflarv$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    yflarv$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(yflarv$CRUISE[i]),conditionMessage(e), "\n")})
}

##PLAICE: 
apeggraw<-read.csv(file='BS_PlaiceEggCatch.csv',header=TRUE,check.names=TRUE)
aplarvraw<-read.csv(file='BS_PlaiceLarvaeCatch.csv',header=TRUE,check.names=TRUE)

apeggraw<-apeggraw[apeggraw$HAUL_ID!='1SS02 81 1 60BON 2',]
aplarvraw<-aplarvraw[aplarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

apeggraw$doy<-as.numeric(mdy.date(apeggraw$MONTH_,apeggraw$DAY_,1960))
aplarvraw$doy<-as.numeric(mdy.date(aplarvraw$MONTH_,aplarvraw$DAY_,1960))
apeggraw$id<-paste(apeggraw$CRUISE,apeggraw$LAT,apeggraw$LON,apeggraw$GMT_DATE_TIME,apeggraw$MESH,sep="_")
aplarvraw$id<-paste(aplarvraw$CRUISE,aplarvraw$LAT,aplarvraw$LON,aplarvraw$GMT_DATE_TIME,aplarvraw$MESH,sep="_")

apeggraw[apeggraw$NET==3,]
apegg_1_2<-apeggraw[apeggraw$NET<3,]
tmp<-table(apegg_1_2$id)
id_2<-names(tmp)[tmp==2]
apegg_1<-apegg_1_2[apegg_1_2$NET==1,]
apegg_2<-apegg_1_2[apegg_1_2$id%in%id_2&apegg_1_2$NET==1,]
apegg_1net<-apegg_1_2[!apegg_1_2$id%in%id_2,]
apegg<-rbind(apegg_1net,apegg_2)
dim(apegg)#We should end up with n = 6181-27-84= 6070 data points

aplarvraw[aplarvraw$NET==3,]
aplarv_1_2<-aplarvraw[aplarvraw$NET<3,]
tmp<-table(aplarv_1_2$id)
id_2<-names(tmp)[tmp==2]
aplarv_2<-aplarv_1_2[aplarv_1_2$id%in%id_2&aplarv_1_2$NET==1,]

aplarv_1net<-aplarv_1_2[!aplarv_1_2$id%in%id_2,]
aplarvae<-rbind(aplarv_1net,aplarv_2)
dim(aplarvae)#We should end up with n = 6181-27-84= 6070 data points

apegg<-apegg[apegg$PRIMARY_NET=='Y',]
aplarvae<-aplarvae[aplarvae$PRIMARY_NET=='Y',]

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
aplarvae$SS<-as.numeric(tmp[index])

attach(apegg)
apegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(apegg)

tmp<-table(apegg$YEAR_)
index<-match(apegg$YEAR_,names(tmp))
apegg$SS<-as.numeric(tmp[index])

apsub<-apegg[apegg$MONTH_>2&apegg$MONTH_<7,]
apsub<-apsub[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]
aplarv<-aplarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]

names(apsub)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')
names(aplarv)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                 'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')

apsub$SSB<-NA
aplarv$SSB<-NA

SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
apSSB<-SSBdf[SSBdf$Species=='alaska plaice',]
apSSB<-apSSB[c('Year','SSB')]

for(i in 1:nrow(aplarv)){
  aplarv$SSB[i]<-apSSB$SSB[apSSB$Year==aplarv$year[i]]
}

for(i in 1:nrow(apsub)){
  apsub$SSB[i]<-apSSB$SSB[apSSB$Year==apsub$year[i]]
}

uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")

idx4<-match(apsub$CRUISE,uniquecruise$unique_cruise)
apsub$ctd_file_name<-uniquecruise$file_name[idx4]

apsub<-na.exclude(apsub)
apsub$temperature<-NA  
apsub$salinity<-NA

for (i in 1:nrow(apsub)){
  tryCatch({
    filename<-paste0(apsub$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==apsub$HAUL_NAME[i]&ctd.file$Station==apsub$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    apsub$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    apsub$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(apsub$CRUISE[i]),conditionMessage(e), "\n")})
}

idx5<-match(aplarv$CRUISE,uniquecruise$unique_cruise)
aplarv$ctd_file_name<-uniquecruise$file_name[idx5]

aplarv<-na.exclude(aplarv)
aplarv$temperature<-NA  
aplarv$salinity<-NA

for (i in 1:nrow(aplarv)){
  tryCatch({
    filename<-paste0(aplarv$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==aplarv$HAUL_NAME[i]&ctd.file$Station==aplarv$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    aplarv$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    aplarv$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(aplarv$CRUISE[i]),conditionMessage(e), "\n")})
}

##FLATHEAD: 
fheggraw<-read.csv(file='BS_FlatheadEggCatch.csv',header=TRUE,check.names=TRUE)
fhlarvraw<-read.csv(file='BS_FlatheadLarvaeCatch.csv',header=TRUE,check.names=TRUE)

fheggraw<-fheggraw[fheggraw$HAUL_ID!='1SS02 81 1 60BON 2',]
fhlarvraw<-fhlarvraw[fhlarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

fheggraw$doy<-as.numeric(mdy.date(fheggraw$MONTH_,fheggraw$DAY_,1960))
fhlarvraw$doy<-as.numeric(mdy.date(fhlarvraw$MONTH_,fhlarvraw$DAY_,1960))
fheggraw$id<-paste(fheggraw$CRUISE,fheggraw$LAT,fheggraw$LON,fheggraw$GMT_DATE_TIME,fheggraw$MESH,sep="_")
fhlarvraw$id<-paste(fhlarvraw$CRUISE,fhlarvraw$LAT,fhlarvraw$LON,fhlarvraw$GMT_DATE_TIME,fhlarvraw$MESH,sep="_")

fheggraw[fheggraw$NET==3,]
fhegg_1_2<-fheggraw[fheggraw$NET<3,]
tmp<-table(fhegg_1_2$id)
id_2<-names(tmp)[tmp==2]
fhegg_1<-fhegg_1_2[fhegg_1_2$NET==1,]
fhegg_2<-fhegg_1_2[fhegg_1_2$id%in%id_2&fhegg_1_2$NET==1,]
fhegg_1net<-fhegg_1_2[!fhegg_1_2$id%in%id_2,]
fhegg<-rbind(fhegg_1net,fhegg_2)
dim(fhegg)#We should end up with n = 6181-27-84= 6070 data points

fhlarvraw[fhlarvraw$NET==3,]
fhlarv_1_2<-fhlarvraw[fhlarvraw$NET<3,]
tmp<-table(fhlarv_1_2$id)
id_2<-names(tmp)[tmp==2]
fhlarv_2<-fhlarv_1_2[fhlarv_1_2$id%in%id_2&fhlarv_1_2$NET==1,]
fhlarv_1net<-fhlarv_1_2[!fhlarv_1_2$id%in%id_2,]
fhlarvae<-rbind(fhlarv_1net,fhlarv_2)
dim(fhlarvae)#We should end up with n = 6181-27-84= 6070 data points

fhegg<-fhegg[fhegg$PRIMARY_NET=='Y',]
fhlarvae<-fhlarvae[fhlarvae$PRIMARY_NET=='Y',]

fhegg$dist<-NA;fhlarvae$dist<-NA
for(i in 1:nrow(fhegg)){
  fhegg$dist[i]<-min(distance.function(fhegg$LAT[i],fhegg$LON[i],fhegg$LAT[fhegg$LARVALCATCHPER10M2>0],fhegg$LON[fhegg$LARVALCATCHPER10M2>0]))
  fhlarvae$dist[i]<-min(distance.function(fhlarvae$LAT[i], fhlarvae$LON[i], fhlarvae$LAT[fhlarvae$LARVALCATCHPER10M2>0], fhlarvae$LON[fhlarvae$LARVALCATCHPER10M2>0]))
}

attach(fhlarvae)
fhlarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(fhlarvae)

tmp<-table(fhlarvae$YEAR_)
index<-match(fhlarvae$YEAR_,names(tmp))
fhlarvae$SS<-as.numeric(tmp[index])

attach(fhegg)
fhegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(fhegg)

tmp<-table(fhegg$YEAR_)
index<-match(fhegg$YEAR_,names(tmp))
fhegg$SS<-as.numeric(tmp[index])

fhsub<-fhegg[fhegg$MONTH_>1&fhegg$MONTH_<7,]
fhsub<-fhsub[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]
fhlarv<-fhlarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','LAT','LON','doy','VOLUME_FILTERED','count','SS')]

names(fhsub)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')
names(fhlarv)<-c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2',
                 'LARVALCATCHPER1000M3','year','lat','lon','doy','vol','count','SS')

fhsub$SSB<-NA
fhlarv$SSB<-NA

SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
fhSSB<-SSBdf[SSBdf$Species=='flathead sole',]
fhSSB<-fhSSB[c('Year','SSB')]

for(i in 1:nrow(fhlarv)){
  fhlarv$SSB[i]<-fhSSB$SSB[fhSSB$Year==fhlarv$year[i]]
}

for(i in 1:nrow(fhsub)){
  fhsub$SSB[i]<-fhSSB$SSB[fhSSB$Year==fhsub$year[i]]
}

uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")

idx6<-match(fhsub$CRUISE,uniquecruise$unique_cruise)
fhsub$ctd_file_name<-uniquecruise$file_name[idx6]

fhsub<-na.exclude(fhsub)
fhsub$temperature<-NA  
fhsub$salinity<-NA

for (i in 1:nrow(fhsub)){
  tryCatch({
    filename<-paste0(fhsub$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==fhsub$HAUL_NAME[i]&ctd.file$Station==fhsub$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    fhsub$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    fhsub$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(fhsub$CRUISE[i]),conditionMessage(e), "\n")})
}

idx7<-match(fhlarv$CRUISE,uniquecruise$unique_cruise)
fhlarv$ctd_file_name<-uniquecruise$file_name[idx7]

fhlarv<-na.exclude(fhlarv)
fhlarv$temperature<-NA  
fhlarv$salinity<-NA

for (i in 1:nrow(fhlarv)){
  tryCatch({
    filename<-paste0(fhlarv$CRUISE[i],".csv")
    ctd.file<-read.csv(file=filename,header=TRUE,check.names=TRUE,na.strings="")
    num.cols<-c(1,2,4:9)
    char.cols<-c(3,12,13,15:17)
    int.cols<-c(10,11,14)
    ctd.file[,num.cols]=lapply(ctd.file[,num.cols],as.numeric)
    ctd.file[,char.cols]=lapply(ctd.file[,char.cols],as.character)
    ctd.file[,int.cols]=lapply(ctd.file[,int.cols],as.integer)
    idx1<-ctd.file$Haul==fhlarv$HAUL_NAME[i]&ctd.file$Station==fhlarv$STATION_NAME[i]
    ctdcast<-ctd.file[idx1,]
    fhlarv$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    fhlarv$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
  }, error=function(e){cat(as.character(print(i)),as.character(fhlarv$CRUISE[i]),conditionMessage(e), "\n")})
}
