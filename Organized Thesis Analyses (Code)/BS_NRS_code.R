
# Northern Rock Sole Script For Loading and Cleaning Data ------------------------
## The following script cleans larval Pacific cod data from the Bering Sea. 
## Northern rock sole spawn in waters difficult to access by ecoFOCI cruises, so egg data are not available. 


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
nrslarvraw<-read.csv(file='BS_NorthernRockSole_Larvae_CatchwZeros_PrimaryNPQ.csv')
setwd(gitwd)

nrslarvraw<-nrslarvraw[nrslarvraw$HAUL_ID!='1SS02 81 1 60BON 2',]

nrslarvraw$doy<-as.numeric(mdy.date(nrslarvraw$MONTH,nrslarvraw$DAY,1960))
nrslarvraw$id<-paste(nrslarvraw$CRUISE,nrslarvraw$LAT,nrslarvraw$LON,nrslarvraw$GMT_DATE_TIME,nrslarvraw$MESH,sep="_")

#net cleaning, larvae; this removes duplicate gear to ensure we're looking at proper catch amounts: 
nrslarvraw[nrslarvraw$NET==3,]
nrslarv_1_2<-nrslarvraw[nrslarvraw$NET<3,]
tmp<-table(nrslarv_1_2$id)
id_2<-names(tmp)[tmp==2]
nrslarv_2<-nrslarv_1_2[nrslarv_1_2$id%in%id_2&nrslarv_1_2$NET==1,]

nrslarv_1net<-nrslarv_1_2[!nrslarv_1_2$id%in%id_2,]
nrslarvae<-rbind(nrslarv_1net,nrslarv_2)
dim(nrslarvae)#We should end up with n = 6181-27-84= 6070 data points

#remove stations where primary net =/= to 
nrslarvae<-nrslarvae[nrslarvae$PRIMARY_NET=='Y',]

#add distance to closest positive catch
nrslarvae$dist<-NA
for(i in 1:nrow(nrslarvae)){
  nrslarvae$dist[i]<-min(distance.function(nrslarvae$LAT[i], 
                                          nrslarvae$LON[i], nrslarvae$LAT[nrslarvae$LARVALCATCHPER10M2>0], 
                                          nrslarvae$LON[nrslarvae$LARVALCATCHPER10M2>0]))
}

attach(nrslarvae)
nrslarvae$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(nrslarvae)

tmp<-table(nrslarvae$YEAR)
index<-match(nrslarvae$YEAR,names(tmp))
nrslarvae$SS<-as.numeric(tmp[index])

nrslarvae$DATE<-NULL
nrslarvae$DATE<-paste(nrslarvae$MONTH,nrslarvae$DAY,nrslarvae$YEAR,sep="/") #put date column in proper format

nrslarv<-nrslarvae[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID',
                   'LARVALCATCHPER10M2','LARVALCATCHPER1000M3','YEAR','MONTH','LAT','LON','doy','VOLUME_FILTERED',
                   'BOTTOM_DEPTH','id','count','SS','DATE')]
nrslarv<-subset(nrslarv,BOTTOM_DEPTH>40&BOTTOM_DEPTH<300)
nrslarv<-subset(nrslarv,LAT<64.5)
names(nrslarv)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','Cper10m2',
                 'Cper1000m3','year','month','lat','lon','doy','vol','bottom_depth','id','count','SS','date')

##Now that these data are properly cleaned, we can add in the CTD data
nrslarv.ctd<-nrslarv

allctd$Link_ID<-NA
nrslarv.ctd$Link_ID<-NA

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier to match 2 sets
nrslarv.ctd$Link_ID<-paste(nrslarv.ctd$CRUISE,nrslarv.ctd$STATION,nrslarv.ctd$HAUL,sep="_")

#now we match CTD data to larval data: 
nrslarv.ctd$temperature<-NA
nrslarv.ctd$salinity<-NA
nrslarv.ctd$CTD_date<-NA
nrslarv.ctd$CTD_link<-NA
nrslarv.ctd$CTD_time<-NA

for(i in 1:nrow(nrslarv.ctd)){
  tryCatch({
    idx1.lv<-allctd$Link_ID==nrslarv.ctd$Link_ID[i]
    ctdcast<-allctd[idx1.lv,]
    nrslarv.ctd$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    nrslarv.ctd$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    nrslarv.ctd$CTD_date[i]<-unique(ctdcast$Date)
    nrslarv.ctd$CTD_link[i]<-unique(ctdcast$Link_ID)
    nrslarv.ctd$CTD_time[i]<-unique(ctdcast$Time)
  },error=function(e){cat(as.character(print(i)),as.character(nrslarv.ctd$CRUISE[i]),
                          conditionMessage(e),"\n")})}
nrslarv.ctd<-nrslarv.ctd[!nrslarv.ctd$temperature=="NaN",]

nrslarv.ctd$date_diff<-NA #now we want to make sure that larval samples were collected 
#at a similar time to when the CTD samples were collected
nrslarv.ctd$date<-parse_date_time(nrslarv.ctd$date,orders="mdy")
nrslarv.ctd$CTD_date<-parse_date_time(nrslarv.ctd$CTD_date,orders="mdy")
nrslarv.ctd$date_diff<-difftime(nrslarv.ctd$date,nrslarv.ctd$CTD_date,units="hour")
nrslarv.ctd<-nrslarv.ctd[which(nrslarv.ctd$date_diff>(-6)&nrslarv.ctd$date_diff<6),]
dim(nrslarv.ctd) #2626, 25

write.csv(nrslarv,'./Ichthyo Data/Cleaned_Cut_NrsLarv.csv')
write.csv(nrslarv.ctd,'./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv')

# Loading in cleaned and subset data, avoiding all above gymnastics -------

nrslarv.ctd<-read.csv('./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv')

nrslarv<-read.csv('./Ichthyo Data/Cleaned_Cut_NrsLarv.csv')

