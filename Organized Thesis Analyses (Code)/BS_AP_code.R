###Alaska Plaice BS Code: 
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

##PLAICE: 

###loading data, subsetting and cleaning properly 
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
table(apegg_1_2$NET)
tmp<-table(apegg_1_2$id)
id_2<-names(tmp)[tmp==2]
apegg_1<-apegg_1_2[apegg_1_2$NET==1,]
apegg_2<-apegg_1_2[apegg_1_2$id%in%id_2&apegg_1_2$NET==1,]

apegg_1net<-apegg_1_2[!apegg_1_2$id%in%id_2,]
apegg<-rbind(apegg_1net,apegg_2)
dim(apegg)#We should end up with n = 6181-27-84= 6070 data points

#net cleaning, larvae: 
aplarvraw[aplarvraw$NET==3,]
aplarv_1_2<-aplarvraw[aplarvraw$NET<3,]
tmp<-table(aplarv_1_2$id)
id_2<-names(tmp)[tmp==2]
aplarv_2<-aplarv_1_2[aplarv_1_2$id%in%id_2&aplarv_1_2$NET==1,]

aplarv_1net<-aplarv_1_2[!aplarv_1_2$id%in%id_2,]
aplarvae<-rbind(aplarv_1net,aplarv_2)
dim(aplarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to Y
apegg<-apegg[apegg$PRIMARY_NET=='Y',]
aplarvae<-aplarvae[aplarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
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

apegg$DATE<-NULL
apegg$DATE<-paste(apegg$MONTH_,apegg$DAY_,apegg$YEAR_,sep="/")

aplarvae$DATE<-NULL
aplarvae$DATE<-paste(aplarvae$MONTH_,aplarvae$DAY_,aplarvae$YEAR_,sep="/")

apsub<-apegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')]
apsub<-subset(apsub,doy>99&doy<182)
aplarv<-aplarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]

names(apsub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(aplarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

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

apsub.ctd<-apsub
aplarv.ctd<-aplarv

#uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")
#idx.eg<-match(apsub.ctd$CRUISE,uniquecruise$unique_cruise)
#idx.lv<-match(aplarv.ctd$CRUISE,uniquecruise$unique_cruise)
#apsub.ctd$ctd_file_name<-uniquecruise$file_name[idx.eg]
#aplarv.ctd$ctd_file_name<-uniquecruise$file_name[idx.lv]

allctd$Link_ID<-NA
apsub.ctd$Link_ID<-NA
aplarv.ctd$Link_ID<-NA

#allctd$Station<-as.numeric(allctd$Station)#conversion in an attempt to make matching work 
#allctd$Haul<-as.numeric(allctd$Haul)
#apsub.ctd$STATION_NAME<-as.numeric(apsub.ctd$STATION) #warnings introduced for cruises: 2MF92, 10C88, 1MP91, 0MF91, 3MF79, 1DN88, MF862 
#apsub.ctd$HAUL_NAME<-as.numeric(apsub.ctd$HAUL)

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
apsub.ctd$Link_ID<-paste(apsub.ctd$CRUISE,apsub.ctd$STATION,apsub.ctd$HAUL,sep="_")
aplarv.ctd$Link_ID<-paste(aplarv.ctd$CRUISE,aplarv.ctd$STATION,aplarv.ctd$HAUL,sep="_")

#apsub.ctd<-apsub.ctd[!is.na(apsub.ctd$ctd_file_name),] 
apsub.ctd$temperature<-NA
apsub.ctd$salinity<-NA
apsub.ctd$CTD_date<-NA
apsub.ctd$CTD_link<-NA

for(i in 1:nrow(apsub.ctd)){
  tryCatch({
    idx1.eg<-allctd$Link_ID==apsub.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.eg,]
    apsub.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    apsub.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    apsub.ctd$CTD_date[i]<-as.character(ctdcast$Date)
    apsub.ctd$CTD_link[i]<-as.character(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(apsub.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
apsub.ctd<-apsub.ctd[!apsub.ctd$temperature=="NaN",] #dim: 1640 24

#aplarv.ctd<-aplarv.ctd[!is.na(aplarv.ctd$ctd_file_name),] 
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
aplarv.ctd<-aplarv.ctd[!aplarv.ctd$temperature=="NaN",]#dim: 2646 24

aplarv.ctd$date_diff<-NA
aplarv.ctd$date<-parse_date_time(aplarv.ctd$date,orders="mdy")
aplarv.ctd$CTD_date<-parse_date_time(aplarv.ctd$CTD_date,orders="mdy")
aplarv.ctd$date_diff<-difftime(aplarv.ctd$date,aplarv.ctd$CTD_date,units="hour")
aplarv.ctd<-aplarv.ctd[which(aplarv.ctd$date_diff>(-6)&aplarv.ctd$date_diff<6),]
dim(aplarv.ctd) #2565

#Additional trimming based on biological characteristics: 
#Trim eggs to a depth less than 150m (spawn on the middle shelf between 50-100m)
#Trim larvae to a depth less than 150m (larvae transported to coastal areas after hatching)
apsub<-subset(apsub,bottom_depth<151)
apsub.ctd<-subset(apsub.ctd,bottom_depth<151)
aplarv<-subset(aplarv,bottom_depth<151)
aplarv.ctd<-subset(aplarv.ctd,bottom_depth<151)

write.csv(apsub,'Cleaned_Cut_ApEggs.csv')
write.csv(apsub.ctd,'Cleaned_Cut_ApEggs_wCTD.csv')
write.csv(aplarv,'Cleaned_Cut_ApLarv.csv')
write.csv(aplarv.ctd,'Cleaned_Cut_ApLarv_wCTD.csv')


###to avoid doing all the above gymnastics: 
apsub<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApEggs.csv',header=TRUE,check.names=TRUE)
apsub.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApEggs_wCTD.csv',header=TRUE,check.names=TRUE)
aplarv<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv.csv',header=TRUE,check.names=TRUE)
aplarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE)

