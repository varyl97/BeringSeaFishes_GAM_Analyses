##Loading NCEP Data for BS: 
#9/2/2021
library(raster)
library(tidyr)

ncep<-brick("C:/Users/varyl/Desktop/MS Research/RStudio/sst.mon.mean 2.nc")
ncep

#plot(ncep[[1]]) #check it out 
#map("world",fill=T,col="wheat2",add=T)

sst.df<-raster::as.data.frame(ncep,xy=TRUE)

sst.df<-subset(sst.df,y>50&y<70)
sst.df<-subset(sst.df,x>179&x<220)
sst.df$x<-sst.df$x-360
sst.df<-na.exclude(sst.df)

ggplot()+geom_tile(data=sst.df,aes(x=x,y=y,color=X1891.01.01,fill=X1891.01.01))+
  geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))

sst.df<-sst.df[,c(1,2,1011:1558)]

feb.sst<-sst.df[,c(1,2,4,16,28,40,52,64,76,88,100,112,124,136,148,160,172,184,196,
                   208,220,232,244,256,268,280,292,304,316,328,340,352,364,
                   376,388,400,412,424,436,448,460,472,484,496,508,520,532,544)]
View(feb.sst)
names(feb.sst)<-c('lon','lat','1975','1976','1977','1978','1979','1980','1981','1982','1983','1984','1985',
                  '1986','1987','1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007','2008',
                  '2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020')

write.csv(feb.sst,"Feb_SST_NCEP_BS.csv")

#Now: 
feb.sst<-read.csv('Feb_SST_NCEP_BS.csv',header=TRUE,check.names=TRUE)
names(feb.sst)<-c('X','lon','lat','1975','1976','1977','1978','1979','1980','1981',
                  '1982','1983','1984','1985','1986','1987',
                  '1988','1989','1990','1991','1992','1993','1994','1995',
                  '1996','1997','1998','1999','2000','2001','2002','2003','2004',
                  '2005','2006','2007','2008','2009','2010','2011','2012','2013',
                  '2014','2015','2016','2017','2018','2019','2020')
feb.sst<-subset(feb.sst,lon<(-157))

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

write.csv(feb.sst,'Feb_SST_ByLocation_NCEP_BS.csv') 
#this dataframe has SST values for a given
#  location averaged for the month of February 

#to create the data frame of regionally averaged Feb SST for each year: 
reg.sst<-aggregate(SST~year,data=feb.sst,FUN=function(feb.sst) (mean=mean(feb.sst)))

write.csv(reg.sst,'Feb_SST_RegionalIndex_NCEP_BS.csv')


#erroneous things: 
sst.df$x<-rep(1:360,180)
sst.df$y<-rev(rep(-89.5:89.5,each=360))
sst.df$cell_id<-1:64800

sst.df<-subset(sst.df,x>179&x<210)
sst.df<-subset(sst.df,y>50&y<75)









