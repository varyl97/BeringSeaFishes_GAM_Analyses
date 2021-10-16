###Alaska Plaice BS Code: 
###CTD Loading: ################################################################
allctd<-read.csv(file="../Environmental Data/All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference

allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,]

##FLATHEAD SOLE: 

###loading data, subsetting and cleaning properly 
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
table(fhegg_1_2$NET)
tmp<-table(fhegg_1_2$id)
id_2<-names(tmp)[tmp==2]
fhegg_1<-fhegg_1_2[fhegg_1_2$NET==1,]
fhegg_2<-fhegg_1_2[fhegg_1_2$id%in%id_2&fhegg_1_2$NET==1,]

fhegg_1net<-fhegg_1_2[!fhegg_1_2$id%in%id_2,]
fhegg<-rbind(fhegg_1net,fhegg_2)
dim(fhegg)#We should end up with n = 6181-27-84= 6070 data points

#net cleaning, larvae: 
fhlarvraw[fhlarvraw$NET==3,]
fhlarv_1_2<-fhlarvraw[fhlarvraw$NET<3,]
tmp<-table(fhlarv_1_2$id)
id_2<-names(tmp)[tmp==2]
fhlarv_2<-fhlarv_1_2[fhlarv_1_2$id%in%id_2&fhlarv_1_2$NET==1,]

fhlarv_1net<-fhlarv_1_2[!fhlarv_1_2$id%in%id_2,]
fhlarvae<-rbind(fhlarv_1net,fhlarv_2)
dim(fhlarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to Y
fhegg<-fhegg[fhegg$PRIMARY_NET=='Y',]
fhlarvae<-fhlarvae[fhlarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
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

fhegg$DATE<-NULL
fhegg$DATE<-paste(fhegg$MONTH_,fhegg$DAY_,fhegg$YEAR_,sep="/")

fhlarvae$DATE<-NULL
fhlarvae$DATE<-paste(fhlarvae$MONTH_,fhlarvae$DAY_,fhlarvae$YEAR_,sep="/")

fhsub<-fhegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')]
fhsub<-subset(fhsub,MONTH_>1&MONTH_<7)
fhlarv<-fhlarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]

names(fhsub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
names(fhlarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

fhsub$SSB<-NA
fhlarv$SSB<-NA

SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
fhSSB<-SSBdf[SSBdf$Species=='alaska plaice',]
fhSSB<-fhSSB[c('Year','SSB')]

for(i in 1:nrow(fhlarv)){
  fhlarv$SSB[i]<-fhSSB$SSB[fhSSB$Year==fhlarv$year[i]]
}

for(i in 1:nrow(fhsub)){
  fhsub$SSB[i]<-fhSSB$SSB[fhSSB$Year==fhsub$year[i]]
}

fhsub.ctd<-fhsub
fhlarv.ctd<-fhlarv

#uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")
#idx.eg<-match(fhsub.ctd$CRUISE,uniquecruise$unique_cruise)
#idx.lv<-match(fhlarv.ctd$CRUISE,uniquecruise$unique_cruise)
#fhsub.ctd$ctd_file_name<-uniquecruise$file_name[idx.eg]
#fhlarv.ctd$ctd_file_name<-uniquecruise$file_name[idx.lv]

allctd$Link_ID<-NA
fhsub.ctd$Link_ID<-NA
fhlarv.ctd$Link_ID<-NA

#allctd$Station<-as.numeric(allctd$Station)#conversion in an attempt to make matching work 
#allctd$Haul<-as.numeric(allctd$Haul)
#fhsub.ctd$STATION_NAME<-as.numeric(fhsub.ctd$STATION) #warnings introduced for cruises: 2MF92, 10C88, 1MP91, 0MF91, 3MF79, 1DN88, MF862 
#fhsub.ctd$HAUL_NAME<-as.numeric(fhsub.ctd$HAUL)

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
fhsub.ctd$Link_ID<-paste(fhsub.ctd$CRUISE,fhsub.ctd$STATION,fhsub.ctd$HAUL,sep="_")
fhlarv.ctd$Link_ID<-paste(fhlarv.ctd$CRUISE,fhlarv.ctd$STATION,fhlarv.ctd$HAUL,sep="_")

#fhsub.ctd<-fhsub.ctd[!is.na(fhsub.ctd$ctd_file_name),] 
fhsub.ctd$temperature<-NA
fhsub.ctd$salinity<-NA
fhsub.ctd$CTD_date<-NA
fhsub.ctd$CTD_link<-NA

for(i in 1:nrow(fhsub.ctd)){
  tryCatch({
    idx1.eg<-allctd$Link_ID==fhsub.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.eg,]
    fhsub.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    fhsub.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    fhsub.ctd$CTD_date[i]<-as.character(ctdcast$Date)
    fhsub.ctd$CTD_link[i]<-as.character(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(fhsub.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
fhsub.ctd<-fhsub.ctd[!fhsub.ctd$temperature=="NaN",] #dim: 1640 24

#fhlarv.ctd<-fhlarv.ctd[!is.na(fhlarv.ctd$ctd_file_name),] 
fhlarv.ctd$temperature<-NA
fhlarv.ctd$salinity<-NA
fhlarv.ctd$CTD_date<-NA
fhlarv.ctd$CTD_link<-NA
fhlarv.ctd$CTD_time<-NA

for(i in 1:nrow(fhlarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==fhlarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    fhlarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    fhlarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    fhlarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    fhlarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    fhlarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(fhlarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
fhlarv.ctd<-fhlarv.ctd[!fhlarv.ctd$temperature=="NaN",]#dim: 2646 24

fhlarv.ctd$date_diff<-NA
fhlarv.ctd$date<-parse_date_time(fhlarv.ctd$date,orders="mdy")
fhlarv.ctd$CTD_date<-parse_date_time(fhlarv.ctd$CTD_date,orders="mdy")
fhlarv.ctd$date_diff<-difftime(fhlarv.ctd$date,fhlarv.ctd$CTD_date,units="hour")
fhlarv.ctd<-fhlarv.ctd[which(fhlarv.ctd$date_diff>(-6)&fhlarv.ctd$date_diff<6),]
dim(fhlarv.ctd) #2565

write.csv(fhsub,'Cleaned_Cut_FhEggs.csv')
write.csv(fhsub.ctd,'Cleaned_Cut_FhEggs_wCTD.csv')
write.csv(fhlarv,'Cleaned_Cut_FhLarv.csv')
write.csv(fhlarv.ctd,'Cleaned_Cut_FhLarv_wCTD.csv')


###to avoid doing all the above gymnastics: 
fhsub<-read.csv(file='../Ichthyo Data/Cleaned_Cut_FhEggs.csv',header=TRUE,check.names=TRUE)
fhsub.ctd<-read.csv(file='../Ichthyo Data/Cleaned_Cut_FhEggs_wCTD.csv',header=TRUE,check.names=TRUE)
fhlarv<-read.csv(file='../Ichthyo Data/Cleaned_Cut_FhLarv.csv',header=TRUE,check.names=TRUE)
fhlarv.ctd<-read.csv(file='../Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE)

