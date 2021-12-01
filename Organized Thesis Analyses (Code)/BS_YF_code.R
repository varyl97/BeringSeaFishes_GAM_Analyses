##Bering Sea Yellowfin Sole: 
#the following code loads and cleans the yellowfin sole egg and larval data from ecoFOCI cruises. 
#at the end of this script are finalized datasets to load in the future. 
###CTD Loading: ################################################################
allctd<-read.csv(file="All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference

allctd<-allctd[allctd$Temperature<14,] #these 2 lines remove anomalous temp/sal values
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,] #this loads a compilation of all CTD data from ecoFOCI trawls 
  #ultimately, these data will be matched with larval data for larval biogeography GAMs. 

###loading data, subsetting and cleaning properly 
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
table(yfegg_1_2$NET)
tmp<-table(yfegg_1_2$id)
id_2<-names(tmp)[tmp==2]
yfegg_1<-yfegg_1_2[yfegg_1_2$NET==1,]
yfegg_2<-yfegg_1_2[yfegg_1_2$id%in%id_2&yfegg_1_2$NET==1,]

yfegg_1net<-yfegg_1_2[!yfegg_1_2$id%in%id_2,]
yfegg<-rbind(yfegg_1net,yfegg_2)
dim(yfegg)#We should end up with n = 6181-27-84= 6070 data points

#net cleaning, larvae: 
yflarvraw[yflarvraw$NET==3,]
yflarv_1_2<-yflarvraw[yflarvraw$NET<3,]
tmp<-table(yflarv_1_2$id)
id_2<-names(tmp)[tmp==2]
yflarv_2<-yflarv_1_2[yflarv_1_2$id%in%id_2&yflarv_1_2$NET==1,]

yflarv_1net<-yflarv_1_2[!yflarv_1_2$id%in%id_2,]
yflarvae<-rbind(yflarv_1net,yflarv_2)
dim(yflarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to Y
yfegg<-yfegg[yfegg$PRIMARY_NET=='Y',]
yflarvae<-yflarvae[yflarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
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

yfegg$DATE<-NULL
yfegg$DATE<-paste(yfegg$MONTH_,yfegg$DAY_,yfegg$YEAR_,sep="/")

yflarvae$DATE<-NULL
yflarvae$DATE<-paste(yflarvae$MONTH_,yflarvae$DAY_,yflarvae$YEAR_,sep="/")

yfsub<-yfegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')]
yfsub<-subset(yfsub,MONTH_>3&MONTH_<10)
yflarv<-yflarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]

names(yfsub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(yflarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

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

yfsub.ctd<-yfsub
yflarv.ctd<-yflarv

#uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")
#idx.eg<-match(yfsub.ctd$CRUISE,uniquecruise$unique_cruise)
#idx.lv<-match(yflarv.ctd$CRUISE,uniquecruise$unique_cruise)
#yfsub.ctd$ctd_file_name<-uniquecruise$file_name[idx.eg]
#yflarv.ctd$ctd_file_name<-uniquecruise$file_name[idx.lv]

allctd$Link_ID<-NA
yfsub.ctd$Link_ID<-NA
yflarv.ctd$Link_ID<-NA

#allctd$Station<-as.numeric(allctd$Station)#conversion in an attempt to make matching work 
#allctd$Haul<-as.numeric(allctd$Haul)
#yfsub.ctd$STATION_NAME<-as.numeric(yfsub.ctd$STATION) #warnings introduced for cruises: 2MF92, 10C88, 1MP91, 0MF91, 3MF79, 1DN88, MF862 
#yfsub.ctd$HAUL_NAME<-as.numeric(yfsub.ctd$HAUL)

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
yfsub.ctd$Link_ID<-paste(yfsub.ctd$CRUISE,yfsub.ctd$STATION,yfsub.ctd$HAUL,sep="_")
yflarv.ctd$Link_ID<-paste(yflarv.ctd$CRUISE,yflarv.ctd$STATION,yflarv.ctd$HAUL,sep="_")

#yfsub.ctd<-yfsub.ctd[!is.na(yfsub.ctd$ctd_file_name),] 
yfsub.ctd$temperature<-NA
yfsub.ctd$salinity<-NA
yfsub.ctd$CTD_date<-NA
yfsub.ctd$CTD_link<-NA

for(i in 1:nrow(yfsub.ctd)){
  tryCatch({
    idx1.eg<-allctd$Link_ID==yfsub.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.eg,]
    yfsub.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    yfsub.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    yfsub.ctd$CTD_date[i]<-as.character(ctdcast$Date)
    yfsub.ctd$CTD_link[i]<-as.character(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(yfsub.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
yfsub.ctd<-yfsub.ctd[!yfsub.ctd$temperature=="NaN",] #dim: 1640 24

#yflarv.ctd<-yflarv.ctd[!is.na(yflarv.ctd$ctd_file_name),] 
yflarv.ctd$temperature<-NA
yflarv.ctd$salinity<-NA
yflarv.ctd$CTD_date<-NA
yflarv.ctd$CTD_link<-NA
yflarv.ctd$CTD_time<-NA

for(i in 1:nrow(yflarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==yflarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    yflarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    yflarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    yflarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    yflarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    yflarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(yflarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
yflarv.ctd<-yflarv.ctd[!yflarv.ctd$temperature=="NaN",]#dim: 2646 24

yflarv.ctd$date_diff<-NA
yflarv.ctd$date<-parse_date_time(yflarv.ctd$date,orders="mdy")
yflarv.ctd$CTD_date<-parse_date_time(yflarv.ctd$CTD_date,orders="mdy")
yflarv.ctd$date_diff<-difftime(yflarv.ctd$date,yflarv.ctd$CTD_date,units="hour")
yflarv.ctd<-yflarv.ctd[which(yflarv.ctd$date_diff>(-6)&yflarv.ctd$date_diff<6),]
dim(yflarv.ctd) #2565

#Further trim based on depth and latitude, based on biology
#Yellowfin sole are coastal spawners, so eggs are restricted to less than 175m
#Larvae are then transported inshore to nursery areas, so larvae restricted to less than 175m
#YFS eggs and larvae are rarely observed north of 61dN, so restricted data to latitudes below 
yfsub<-subset(yfsub,bottom_depth<176&lat<62)
yfsub.ctd<-subset(yfsub.ctd,bottom_depth<176&lat<62)
yflarv<-subset(yflarv,bottom_depth<176)
yflarv.ctd<-subset(yflarv.ctd,bottom_depth<176)

write.csv(yfsub,'./Ichthyo Data/Cleaned_Cut_YfEggs.csv')
write.csv(yfsub.ctd,'./Ichthyo Data/Cleaned_Cut_YfEggs_wCTD.csv')
write.csv(yflarv,'./Ichthyo Data/Cleaned_Cut_YfLarv.csv')
write.csv(yflarv.ctd,'./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv')


###to avoid doing all the above gymnastics: 
yfsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfEggs.csv',header=TRUE,check.names=TRUE)
yfsub.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfEggs_wCTD.csv',header=TRUE,check.names=TRUE)
yflarv<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv.csv',header=TRUE,check.names=TRUE)
yflarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv',header=TRUE,check.names=TRUE)



