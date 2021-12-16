##Rex Sole: 
#the following code loads and cleans the rex sole egg and larval data from ecoFOCI cruises. 
#at the end of this script are finalized datasets to load in the future. 


# CTD Loading:  -----------------------------------------------------------
#this loads a compilation of all CTD data from ecoFOCI cruises. 
#these data will be matched with larval data for larval biogeography analyses. 
allctd<-read.csv(file="./Environmental Data/All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference

allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>29&allctd$Salinity<36,] #this removes anomalous temp/sal values


# Loading Data, Subsetting, and Cleaning ----------------------------------
#the following lines related to cleaning are based off intel from collaborators familiar with cruises
rxeggraw<-read.csv(file='BS_RexSole_Eggs_CatchwZeros.csv')
#rxlarvraw<-read.csv(file='BS_RexSole_Larvae_CatchwZeros_PrimaryNPQ.csv')

rxeggraw<-rxeggraw[rxeggraw$HAUL_ID!='1SS02 81 1 60BON 2',]
#rxlarvraw<-rxlarvraw[rxlarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

rxeggraw$doy<-as.numeric(mdy.date(rxeggraw$MONTH,rxeggraw$DAY,1960))
#rxlarvraw$doy<-as.numeric(mdy.date(rxlarvraw$MONTH,rxlarvraw$DAY,1960))
rxeggraw$id<-paste(rxeggraw$CRUISE,rxeggraw$LAT,rxeggraw$LON,rxeggraw$GMT_DATE_TIME,rxeggraw$MESH,sep="_")
#rxlarvraw$id<-paste(rxlarvraw$CRUISE,rxlarvraw$LAT,rxlarvraw$LON,rxlarvraw$GMT_DATE_TIME,rxlarvraw$MESH,sep="_")

rxeggraw[rxeggraw$NET==3,] 
rxegg_1_2<-rxeggraw[rxeggraw$NET<3,]
table(rxegg_1_2$NET)
tmp<-table(rxegg_1_2$id)
id_2<-names(tmp)[tmp==2]
rxegg_1<-rxegg_1_2[rxegg_1_2$NET==1,]
rxegg_2<-rxegg_1_2[rxegg_1_2$id%in%id_2&rxegg_1_2$NET==1,] #this removes duplicate gear, ensuring we are looking at the data
                                                    #with the proper gear type 

rxegg_1net<-rxegg_1_2[!rxegg_1_2$id%in%id_2,]
rxegg<-rbind(rxegg_1net,rxegg_2)
dim(rxegg)#slightly different dimensions than other species, will end up with 6110 rows


#net cleaning, larvae: 
#rxlarvraw[rxlarvraw$NET==3,]
#rxlarv_1_2<-rxlarvraw[rxlarvraw$NET<3,]
#tmp<-table(rxlarv_1_2$id)
#id_2<-names(tmp)[tmp==2]
#rxlarv_2<-rxlarv_1_2[rxlarv_1_2$id%in%id_2&rxlarv_1_2$NET==1,]

#rxlarv_1net<-rxlarv_1_2[!rxlarv_1_2$id%in%id_2,]
#rxlarvae<-rbind(rxlarv_1net,rxlarv_2)
#dim(rxlarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to Y
rxegg<-rxegg[rxegg$PRIMARY_NET=='Y',]
#rxlarvae<-rxlarvae[rxlarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
rxegg$dist<-NA;#rxlarvae$dist<-NA
for(i in 1:nrow(rxegg)){
  rxegg$dist[i]<-min(distance.function(rxegg$LAT[i],rxegg$LON[i],rxegg$LAT[rxegg$LARVALCATCHPER10M2>0],rxegg$LON[rxegg$LARVALCATCHPER10M2>0]))
  #rxlarvae$dist[i]<-min(distance.function(rxlarvae$LAT[i], rxlarvae$LON[i], rxlarvae$LAT[rxlarvae$LARVALCATCHPER10M2>0], rxlarvae$LON[rxlarvae$LARVALCATCHPER10M2>0]))
}

#attach(rxlarvae)
#rxlarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
#detach(rxlarvae)

#tmp<-table(rxlarvae$YEAR)
#index<-match(rxlarvae$YEAR,names(tmp))
#rxlarvae$SS<-as.numeric(tmp[index])

attach(rxegg)
rxegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(rxegg)

tmp<-table(rxegg$YEAR)
index<-match(rxegg$YEAR,names(tmp))
rxegg$SS<-as.numeric(tmp[index])

rxegg$DATE<-NULL
rxegg$DATE<-paste(rxegg$MONTH,rxegg$DAY,rxegg$YEAR,sep="/")

#rxlarvae$DATE<-NULL
#rxlarvae$DATE<-paste(rxlarvae$MONTH,rxlarvae$DAY,rxlarvae$YEAR,sep="/")

rxsub<-rxegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
               'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR','MONTH','LAT','LON','doy','VOLUME_FILTERED',
               'BOTTOM_DEPTH','id','count','SS','DATE')]
rxsub<-subset(rxsub,MONTH>3&MONTH<8)
#rxlarvae<-subset(rxlarvae,LAT<62.1)
#rxlarv<-rxlarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   #'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR','MONTH','LAT','LON','doy','VOLUME_FILTERED',
                   #'BOTTOM_DEPTH','id','count','SS','DATE')]

names(rxsub)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')
#names(rxlarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 #'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

##Now that these data are properly cleaned, we can add in the CTD data
rxsub.ctd<-rxsub
#rxlarv.ctd<-rxlarv

allctd$Link_ID<-NA
rxsub.ctd$Link_ID<-NA
#rxlarv.ctd$Link_ID<-NA

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier to match
rxsub.ctd$Link_ID<-paste(rxsub.ctd$CRUISE,rxsub.ctd$STATION,rxsub.ctd$HAUL,sep="_")
#rxlarv.ctd$Link_ID<-paste(rxlarv.ctd$CRUISE,rxlarv.ctd$STATION,rxlarv.ctd$HAUL,sep="_")

#add CTD to egg data first
rxsub.ctd$temperature<-NA
rxsub.ctd$salinity<-NA
rxsub.ctd$CTD_date<-NA
rxsub.ctd$CTD_link<-NA

for(i in 1:nrow(rxsub.ctd)){
  tryCatch({
    idx1.eg<-allctd$Link_ID==rxsub.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.eg,]
    rxsub.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    rxsub.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    rxsub.ctd$CTD_date[i]<-as.character(ctdcast$Date)
    rxsub.ctd$CTD_link[i]<-as.character(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(rxsub.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
rxsub.ctd<-rxsub.ctd[!rxsub.ctd$temperature=="NaN",] #dim: 1699 23

#rxlarv.ctd$temperature<-NA
#rxlarv.ctd$salinity<-NA
#rxlarv.ctd$CTD_date<-NA
#rxlarv.ctd$CTD_link<-NA
#rxlarv.ctd$CTD_time<-NA

#for(i in 1:nrow(rxlarv.ctd)){
 # tryCatch({
  #  idx1.lv<-allctd$Link_ID==rxlarv.ctd$Link_ID[i]
   # ctdcast<-allctd[idx1.lv,]
    #rxlarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    #rxlarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    #rxlarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    #rxlarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    #rxlarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  #},error=function(e){cat(as.character(print(i)),as.character(rxlarv.ctd$CRUISE[i]),
   #                       conditionMessage(e),"\n")})}
#rxlarv.ctd<-rxlarv.ctd[!rxlarv.ctd$temperature=="NaN",]

#rxlarv.ctd$date_diff<-NA #now we want to make sure that larval samples were collected 
                            #at a similar time to when the CTD samples were collected
#rxlarv.ctd$date<-parse_date_time(rxlarv.ctd$date,orders="mdy")
#rxlarv.ctd$CTD_date<-parse_date_time(rxlarv.ctd$CTD_date,orders="mdy")
#rxlarv.ctd$date_diff<-difftime(rxlarv.ctd$date,rxlarv.ctd$CTD_date,units="hour")
#rxlarv.ctd<-rxlarv.ctd[which(rxlarv.ctd$date_diff>(-6)&rxlarv.ctd$date_diff<6),]
#dim(rxlarv.ctd) 

write.csv(rxsub,'./Ichthyo Data/Cleaned_Cut_RxEggs.csv')
write.csv(rxsub.ctd,'./Ichthyo Data/Cleaned_Cut_RxEggs_wCTD.csv')
#write.csv(rxlarv,'Cleaned_Cut_RxLarv.csv')
#write.csv(rxlarv.ctd,'Cleaned_Cut_RxLarv_wCTD.csv')


# Loading in cleaned and subset data, avoiding all above gymnastics -------

rxsub<-read.csv('./Ichthyo Data/Cleaned_Cut_RxEggs.csv')
#rxlarv.ctd<-read.csv('./Ichthyo Data/Cleaned_Cut_RxLarv_wCTD.csv')

rxsub.ctd<-read.csv('./Ichthyo Data/Cleaned_Cut_RxEggs_wCTD.csv')
#rxlarv<-read.csv('./Ichthyo Data/Cleaned_Cut_RxLarv.csv')






































