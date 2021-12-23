#BS Pollock Code: 

###BATHYMETRY###################################################################
setwd(rwd)
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
setwd(gitwd)
allctd<-read.csv(file="./Environmental Data/All_CTD_Data_8302021.csv")
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
setwd(rwd)
pkeggraw<-read.csv(file='BS_PollockEggCatch.csv',header=TRUE,check.names=TRUE)
pklarvraw<-read.csv(file='BS_PollockLarvaeCatch.csv',header=TRUE,check.names=TRUE)
setwd(gitwd)

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
pksub<-pksub[pksub$LAT>54&pksub$LON>(-173)&pksub$LAT<62,]
pklarv<-pklarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]
pklarv<-pklarv[pklarv$LAT>53.5&pklarv$LON>(-175),]
names(pksub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(pklarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

pksub$SSB<-NA
pklarv$SSB<-NA

SSBdf<-read.csv(file='./Ichthyo Data/SSB_allsp.csv',header=TRUE,check.names=TRUE)
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

#Save for later: 
write.csv(pksub,'./Ichthyo Data/Cleaned_Cut_PkEggs.csv')
write.csv(pksub.ctd,'./Ichthyo Data/Cleaned_Cut_PkEggs_wCTD.csv')
write.csv(pklarv,'./Ichthyo Data/Cleaned_Cut_PkLarv.csv')
write.csv(pklarv.ctd,'./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv')

###to avoid doing all the above gymnastics: if doesn't work, try '../' 
pksub<-read.csv(file='../Ichthyo Data/Cleaned_Cut_PkEggs.csv',header=TRUE,check.names=TRUE)
pksub.ctd<-read.csv(file='../Ichthyo Data/Cleaned_Cut_PkEggs_wCTD.csv',header=TRUE,check.names=TRUE)
pklarv<-read.csv(file='../Ichthyo Data/Cleaned_Cut_PkLarv.csv',header=TRUE,check.names=TRUE)
pklarv.ctd<-read.csv(file='../Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',header=TRUE,check.names=TRUE)
