##Loading NCEP Data for BS: 
#9/2/2021
library(raster)
library(tidyr)

ncep<-brick("./Environmental Data/sst.mon.mean 2.nc")#if this doesn't work, try "../" 
ncep

#plot(ncep[[1]]) #check it out 
#map("world",fill=T,col="wheat2",add=T)

sst.df<-raster::as.data.frame(ncep,xy=TRUE)

sst.df<-subset(sst.df,y>55&y<58)
sst.df<-subset(sst.df,x>(191)&x<(197)) #want -164 to -168
sst.df$x<-sst.df$x-360


sst.df<-na.exclude(sst.df)

world<-map_data("world")
BSmap<-world[world$long<(-155)&world$lat>50&world$lat<70,]
windows()
ggplot()+geom_tile(data=sst.df,aes(x=x,y=y,color=X1891.01.01,fill=X1891.01.01))+
  geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region)) #to check whats going on 

sst.df<-sst.df[,c(1,2,1011:1558)] #trims to years of interest

View(sst.df)#note columns we want to set apart 
feb.sst<-sst.df[,c(1,2,4,16,28,40,52,64,76,88,100,112,124,136,148,160,172,184,196,
                   208,220,232,244,256,268,280,292,304,316,328,340,352,364,
                   376,388,400,412,424,436,448,460,472,484,496,508,520,532,544)]
mar.sst<-sst.df[,c(1,2,5,17,29,41,53,65,77,89,101,113,125,137,149,161,173,185,197,
                   209,221,233,245,257,269,281,293,305,317,329,341,353,365,377,
                   389,401,413,425,437,449,461,473,485,497,509,521,533,545)]#for the spring spawners
may.sst<-sst.df[,c(1,2,7,19,31,43,55,67,79,91,103,115,127,139,151,163,175,187,199,
                   211,223,235,247,259,271,283,295,307,319,331,343,355,367,379,
                   391,403,415,427,439,451,463,475,487,499,511,523,535,547)]#for the summer spawners

View(feb.sst)
names(feb.sst)<-c('lon','lat','1975','1976','1977','1978','1979','1980','1981','1982','1983','1984','1985',
                  '1986','1987','1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008',
                  '2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020')
names(mar.sst)<-c('lon','lat','1975','1976','1977','1978','1979','1980','1981','1982','1983','1984','1985',
                  '1986','1987','1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008',
                  '2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020')
names(may.sst)<-c('lon','lat','1975','1976','1977','1978','1979','1980','1981','1982','1983','1984','1985',
                  '1986','1987','1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008',
                  '2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020')

write.csv(feb.sst,"./Environmental Data/Feb_SST_NCEP_BS.csv")
write.csv(mar.sst,"./Environmental Data/Mar_SST_NCEP_BS.csv")
write.csv(may.sst,"./Environmental Data/May_SST_NCEP_BS.csv")


# Create regional indices with averaged SSTs as well as local indices --------

feb.sst<-read.csv('./Environmental Data/Feb_SST_NCEP_BS.csv',header=TRUE,check.names=TRUE)
names(feb.sst)<-c('X','lon','lat','1975','1976','1977','1978','1979','1980','1981',
                  '1982','1983','1984','1985','1986','1987',
                  '1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004',
                  '2005','2006','2007','2008','2009','2010','2011','2012','2013',
                  '2014','2015','2016','2017','2018','2019','2020')
feb.sst<-subset(feb.sst,lon<(-157)&lat<62&lat>51) #isolate Southeast BS region 

mar.sst<-read.csv('./Environmental Data/Mar_SST_NCEP_BS.csv',header=TRUE,check.names=TRUE)
names(mar.sst)<-c('X','lon','lat','1975','1976','1977','1978','1979','1980','1981',
                  '1982','1983','1984','1985','1986','1987',
                  '1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004',
                  '2005','2006','2007','2008','2009','2010','2011','2012','2013',
                  '2014','2015','2016','2017','2018','2019','2020')
mar.sst<-subset(mar.sst,lon<(-157)&lat<62&lat>51)


may.sst<-read.csv('./Environmental Data/May_SST_NCEP_BS.csv',header=TRUE,check.names=TRUE)
names(may.sst)<-c('X','lon','lat','1975','1976','1977','1978','1979','1980','1981',
                  '1982','1983','1984','1985','1986','1987',
                  '1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004',
                  '2005','2006','2007','2008','2009','2010','2011','2012','2013',
                  '2014','2015','2016','2017','2018','2019','2020')
may.sst<-subset(may.sst,lon<(-157)&lat<62&lat>51)

feb.sst<-gather(feb.sst,year,SST,c('1975','1976','1977','1978','1979','1980',
                                           '1981','1982','1983','1984','1985',
                                           '1986','1987','1988','1989','1990',
                                           '1991','1992','1993','1994','1995',
                                           '1996','1997','1998','1999','2000',
                                           '2001','2002','2003','2004',
                                           '2005','2006','2007','2008','2009',
                                           '2010','2011','2012','2013',
                                           '2014','2015','2016','2017','2018',
                                           '2019','2020'),na.rm=FALSE,
                                            convert=FALSE)
feb.sst<-feb.sst[c('lon','lat','year','SST')]


mar.sst<-gather(mar.sst,year,SST,c('1975','1976','1977','1978','1979','1980',
                                   '1981','1982','1983','1984','1985',
                                   '1986','1987','1988','1989','1990',
                                   '1991','1992','1993','1994','1995',
                                   '1996','1997','1998','1999','2000',
                                   '2001','2002','2003','2004',
                                   '2005','2006','2007','2008','2009',
                                   '2010','2011','2012','2013',
                                   '2014','2015','2016','2017','2018',
                                   '2019','2020'),na.rm=FALSE,
                convert=FALSE)
mar.sst<-mar.sst[c('lon','lat','year','SST')]

may.sst<-gather(may.sst,year,SST,c('1975','1976','1977','1978','1979','1980',
                                   '1981','1982','1983','1984','1985',
                                   '1986','1987','1988','1989','1990',
                                   '1991','1992','1993','1994','1995',
                                   '1996','1997','1998','1999','2000',
                                   '2001','2002','2003','2004',
                                   '2005','2006','2007','2008','2009',
                                   '2010','2011','2012','2013',
                                   '2014','2015','2016','2017','2018',
                                   '2019','2020'),na.rm=FALSE,
                convert=FALSE)
may.sst<-may.sst[c('lon','lat','year','SST')]

write.csv(feb.sst,'./Environmental Data/Feb_SST_ByLocation_NCEP_BS.csv') 
write.csv(mar.sst,'./Environmental Data/Mar_SST_ByLocation_NCEP_BS.csv')
write.csv(may.sst,'./Environmental Data/May_SST_ByLocation_NCEP_BS.csv')
#this dataframe has SST values for a given
#  location averaged for the month of February 

#to create the data frame of regionally averaged Feb SST for each year: 
feb.reg.sst<-aggregate(SST~year,data=feb.sst,FUN=function(feb.sst) (mean=mean(feb.sst)))
write.csv(feb.reg.sst,'./Environmental Data/Feb_SST_RegionalIndex_NCEP_BS.csv')

mar.reg.sst<-aggregate(SST~year,data=mar.sst,FUN=function(mar.sst)(mean=mean(mar.sst)))
write.csv(mar.reg.sst,'./Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv')

may.reg.sst<-aggregate(SST~year,data=may.sst,FUN=function(may.sst)(mean=mean(may.sst)))
write.csv(may.reg.sst,'./Environmental Data/May_SST_RegionalIndex_NCEP_BS.csv')

windows()
ggplot()+geom_tile(data=feb.sst,aes(x=lon,y=lat,fill=SST))+
  scale_fill_viridis_c()+geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=expression(paste("Longitude ("^0,'W)')),
       y=expression(paste("Latitude ("^0,'N)')),)

windows()
ggplot()+geom_tile(data=sub1,aes(x=lon,y=lat,fill=SST))+
  scale_fill_viridis_c()+geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x=expression(paste("Longitude ("^0,'W)')),
       y=expression(paste("Latitude ("^0,'N)')))




