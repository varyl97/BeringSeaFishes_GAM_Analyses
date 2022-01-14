
# Pacific Cod Script For Loading and Cleaning Data ------------------------
## The following script cleans larval Pacific cod data from the Bering Sea. 
## As Pacific cod spawn beneath sea ice, egg data are not available for this species. 


# Load CTD Data -----------------------------------------------------------

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


# Load Larval Data and Clean ----------------------------------------------

#the following lines related to cleaning are based off intel from collaborators familiar with cruises
setwd(rwd)
pclarvraw<-read.csv(file='BS_Pacific_Cod_Larvae_CatchwZeros.csv')
setwd(gitwd)

pclarvraw<-pclarvraw[pclarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

pclarvraw$doy<-as.numeric(mdy.date(pclarvraw$MONTH,pclarvraw$DAY,1960))
pclarvraw$id<-paste(pclarvraw$CRUISE,pclarvraw$LAT,pclarvraw$LON,pclarvraw$GMT_DATE_TIME,pclarvraw$MESH,sep="_")

#net cleaning, larvae; this removes duplicate gear to ensure we're looking at proper catch amounts: 
pclarvraw[pclarvraw$NET==3,]
pclarv_1_2<-pclarvraw[pclarvraw$NET<3,]
tmp<-table(pclarv_1_2$id)
id_2<-names(tmp)[tmp==2]
pclarv_2<-pclarv_1_2[pclarv_1_2$id%in%id_2&pclarv_1_2$NET==1,]

pclarv_1net<-pclarv_1_2[!pclarv_1_2$id%in%id_2,]
pclarvae<-rbind(pclarv_1net,pclarv_2)
dim(pclarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to 
pclarvae<-pclarvae[pclarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
pclarvae$dist<-NA
for(i in 1:nrow(pclarvae)){
  pclarvae$dist[i]<-min(distance.function(pclarvae$LAT[i], 
                                          pclarvae$LON[i], pclarvae$LAT[pclarvae$LARVALCATCHPER10M2>0], 
                                          pclarvae$LON[pclarvae$LARVALCATCHPER10M2>0]))
}

attach(pclarvae)
pclarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(pclarvae)

tmp<-table(pclarvae$YEAR)
index<-match(pclarvae$YEAR,names(tmp))
pclarvae$SS<-as.numeric(tmp[index])

pclarvae$DATE<-NULL
pclarvae$DATE<-paste(pclarvae$MONTH,pclarvae$DAY,pclarvae$YEAR,sep="/") #put date column in proper format

pclarv<-pclarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR','MONTH','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]
pclarv<-subset(pclarv,BOTTOM_DEPTH>40&BOTTOM_DEPTH<250)
names(pclarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

##Now that these data are properly cleaned, we can add in the CTD data
pclarv.ctd<-pclarv

allctd$Link_ID<-NA
pclarv.ctd$Link_ID<-NA

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier to match 2 sets
pclarv.ctd$Link_ID<-paste(pclarv.ctd$CRUISE,pclarv.ctd$STATION,pclarv.ctd$HAUL,sep="_")

#now we match CTD data to larval data: 
pclarv.ctd$temperature<-NA
pclarv.ctd$salinity<-NA
pclarv.ctd$CTD_date<-NA
pclarv.ctd$CTD_link<-NA
pclarv.ctd$CTD_time<-NA

for(i in 1:nrow(pclarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==pclarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    pclarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    pclarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    pclarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    pclarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    pclarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(pclarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
pclarv.ctd<-pclarv.ctd[!pclarv.ctd$temperature=="NaN",]

pclarv.ctd$date_diff<-NA #now we want to make sure that larval samples were collected 
#at a similar time to when the CTD samples were collected
pclarv.ctd$date<-parse_date_time(pclarv.ctd$date,orders="mdy")
pclarv.ctd$CTD_date<-parse_date_time(pclarv.ctd$CTD_date,orders="mdy")
pclarv.ctd$date_diff<-difftime(pclarv.ctd$date,pclarv.ctd$CTD_date,units="hour")
pclarv.ctd<-pclarv.ctd[which(pclarv.ctd$date_diff>(-6)&pclarv.ctd$date_diff<6),]
dim(pclarv.ctd) #2077, 25

pclarv.ctd<-subset(pclarv.ctd,lat<62)

write.csv(pclarv,'./Ichthyo Data/Cleaned_Cut_PcLarv.csv')
write.csv(pclarv.ctd,'./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv')

# Loading in cleaned and subset data, avoiding all above gymnastics -------

pclarv.ctd<-read.csv('./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv')

pclarv<-read.csv('./Ichthyo Data/Cleaned_Cut_PcLarv.csv')

