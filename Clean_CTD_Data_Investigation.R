######################Data Investigation and Checking: 
##Loading in pkegg data to have
pkeggraw<-read.csv(file='BS_PollockEggCatch.csv',header=TRUE,check.names=TRUE)

pkeggraw<-pkeggraw[pkeggraw$HAUL_ID!='1SS02 81 1 60BON 2',]

pkeggraw$doy<-as.numeric(mdy.date(pkeggraw$MONTH_,pkeggraw$DAY_,1960))

pkeggraw$id<-paste(pkeggraw$CRUISE,pkeggraw$LAT,pkeggraw$LON,pkeggraw$GMT_DATE_TIME,pkeggraw$MESH,sep="_")

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

#remove stations where primary net =/= to Y
pkegg<-pkegg[pkegg$PRIMARY_NET=='Y',]

#add distance to closest positive catch
pkegg$dist<-NA
for(i in 1:nrow(pkegg)){
  pkegg$dist[i]<-min(distance.function(pkegg$LAT[i],pkegg$LON[i],pkegg$LAT[pkegg$LARVALCATCHPER10M2>0],pkegg$LON[pkegg$LARVALCATCHPER10M2>0]))
}

attach(pkegg)
pkegg$count<-round(LARVALCATCHPER1000M3*VOLUME_FILTERED/1000,0)
detach(pkegg)

tmp<-table(pkegg$YEAR_)
index<-match(pkegg$YEAR_,names(tmp))
pkegg$SS<-as.numeric(tmp[index])

pkegg$DATE<-NULL
pkegg$DATE<-paste(pkegg$MONTH_,pkegg$DAY_,pkegg$YEAR_,sep="/")

pkegg<-pkegg[c('CRUISE','STATION_NAME','HAUL_NAME','GMT_DATE_TIME','HAUL_ID','LARVALCATCHPER10M2','LARVALCATCHPER1000M3',
               'YEAR_','MONTH_','LAT','LON','doy','VOLUME_FILTERED','BOTTOM_DEPTH','id','count','SS','DATE')]
names(pkegg)<-c('CRUISE','STATION','HAUL','GMT_DATE_TIME','HAUL_ID','cpue10m2','cpue1000m3','year',
                'month','lat','lon','doy','vol','bottom_depth','id','count','SS','DATE')

pkegg$SSB<-NA
SSBdf<-read.csv(file='SSB_allsp.csv',header=TRUE,check.names=TRUE)
SSBdf$SSB<-as.numeric(SSBdf$SSB)
pkSSB<-SSBdf[SSBdf$Species=='walleye pollock',]
pkSSB<-pkSSB[c('Year','SSB')]

for(i in 1:nrow(pkegg)){
  pkegg$SSB[i]<-pkSSB$SSB[pkSSB$Year==pkegg$year[i]]
}

##Combining CTD files into one data frame: 
setwd("C:/Users/varyl/Desktop/MS Research/RStudio/Cruise Working Directory")
cruisewd="C:/Users/varyl/Desktop/MS Research/RStudio/Cruise Working Directory/Cruises (most recent versions)"

list_files = list.files(file.path('Cruises (most recent versions)'))
setwd(cruisewd)
scalar_files = list_files
all_scalar_dfs = lapply(scalar_files, read.csv)

base_scalar=NULL

for(i in 1:length(all_scalar_dfs)){
  filename<-scalar_files[i]
  tmp1<-read.csv(file=filename)
  base_scalar=rbind(base_scalar,tmp1)
}

write.csv(base_scalar,"All_CTD_Data_8302021.csv")

##Now working with the one CTD dataframe/file: 
allctd<-read.csv(file="All_CTD_Data_8302021.csv")
names(allctd)
allctd<-allctd[c('Latitude','Longitude','Date','Time','Pressure','Depth',
                 'Temperature','Conductivity','Salinity','Sigma.T',
                 'Flag','Cruise','Station','Haul','Grid','FOCI_HAUL_ID',
                 'FOCI_file')]

allctd$Year<-NA
allctd$Year<-str_sub(allctd$Date,start=-4) #to get easy year reference

allctd<-allctd[allctd$Temperature<14,]
allctd<-allctd[allctd$Salinity>28|allctd$Salinity<36,]

##Linking CTD to eggs: 
uniquecruise<-read.csv(file='Unique_Cruise_Files.csv',header=TRUE,check.names=TRUE,na.strings="")
idx<-match(pkegg$CRUISE,uniquecruise$unique_cruise)
pkegg$ctd_file_name<-uniquecruise$file_name[idx]

allctd$Link_ID=NULL
pkegg$Link_ID=NULL

allctd$Station<-as.numeric(allctd$Station)#conversion in an attempt to make matching work 
allctd$Haul<-as.numeric(allctd$Haul)
#pkegg$STATION<-as.numeric(pkegg$STATION) #warnings introduced for cruises: 2MF92, 10C88, 1MP91, 0MF91, 3MF79, 1DN88, MF862 
#pkegg$HAUL<-as.numeric(pkegg$HAUL)

allctd$Link_ID<-paste(allctd$Cruise,allctd$Station,allctd$Haul,sep="_") #create new, clean unique identifier
pkegg$Link_ID<-paste(pkegg$CRUISE,pkegg$STATION,pkegg$HAUL,sep="_")

pkegg<-pkegg[!is.na(pkegg$ctd_file_name),] 
pkegg$temperature<-NA
pkegg$salinity<-NA
pkegg$CTD_date<-NA
pkegg$CTD_link<-NA

for(i in 1:nrow(pkegg)){
  tryCatch({
    idx1<-allctd$Link_ID==pkegg$Link_ID[i]
    ctdcast<-allctd[idx1,]
    pkegg$temperature[i]<-mean(ctdcast$Temperature[ctdcast$Depth<11],na.rm=T)
    pkegg$salinity[i]<-mean(ctdcast$Salinity[ctdcast$Depth<11],na.rm=T)
    pkegg$CTD_date[i]<-unique(ctdcast$Date)
    pkegg$CTD_link[i]<-unique(ctdcast$Link_ID)
  },error=function(e){cat(as.character(print(i)),as.character(pkegg$CRUISE[i]),
                          conditionMessage(e),"\n")})}

pkegg$Date_Diff<-NULL
pkegg$DATE<-parse_date_time(pkegg$DATE,orders="mdy") #in lubridate package, best fxn found so far
pkegg$CTD_date<-parse_date_time(pkegg$CTD_date,orders="mdy")
pkegg$Date_Diff<-difftime(pkegg$DATE,pkegg$CTD_date,units="hour")

#this above method isn't quite working because some of the stations/hauls and link 
#         IDs from allctd, despite existing, aren't copying over 
surf.ctd<-data.frame(Latitude=numeric(),Longitude=numeric(),Date=character(),
                     Time=character(),Year=character(),Cruise=character(),
                     Station=numeric(),Haul=numeric(),Top_Temp=numeric(),
                     Top_Sal=numeric(),stringsAsFactors=FALSE)

for(i in 1:nrow(allctd)){
  surf.ctd$Latitude<-as.numeric(allctd$Latitude[i])
  surf.ctd$Longitude[i]<-allctd$Longitude[i]
  surf.ctd$Date[i]<-allctd$Date[i]
  surf.ctd$Time[i]<-allctd$Time[i]
  surf.ctd$Year[i]<-allctd$Year[i]
  surf.ctd$Cruise[i]<-allctd$Cruise[i]
  surf.ctd$Station[i]<-allctd$Station[i]
  surf.ctd$Haul[i]<-allctd$Haul[i]
  surf.ctd$Top_Temp[i]<-mean(allctd$Temperature[allctd$Depth<11],na.rm=T)
  surf.ctd$Top_Sal[i]<-mean(allctd$Salinity[allctd$Depth<11],na.rm=T)
}

#Date checking for larvae: 
pklarv.ctd$date_diff<-NA
pklarv.ctd$date<-parse_date_time(pklarv.ctd$date,orders="mdy")
pklarv.ctd$CTD_date<-parse_date_time(pklarv.ctd$CTD_date,orders="mdy")
pklarv.ctd$date_diff<-difftime(pklarv.ctd$date,pklarv.ctd$CTD_date,units="hour")














