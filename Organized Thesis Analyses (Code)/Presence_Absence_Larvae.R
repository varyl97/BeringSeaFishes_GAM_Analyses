# Generate Presence-Absence Models for Larval Biogeography ----------------
#the following code will generate presence (1) or absence (0) indices for all lines of larval data
#then, a generalized additive model is fit using a binomial distribution 
#this code includes both straight model outputs and predictions based on model (to improve figure visualization)

str_name<-"./Environmental Data/expanded_BS_bathy.tif"
bathy<-raster(str_name) 

# Applying Presence-Absence Values to Data --------------------------------

aplarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE)
fhlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE)
nrslarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pclarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pklarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',header=TRUE,check.names=TRUE)
rxlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_RxLarv_wCTD.csv',header=TRUE,check.names=TRUE)
yflarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv',header=TRUE,check.names=TRUE)

aplarv.ctd <- aplarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                    Cper10m2==0 ~ as.numeric('0')))
fhlarv.ctd <- fhlarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2==0 ~ as.numeric('0')))
nrslarv.ctd <- nrslarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                    Cper10m2==0 ~ as.numeric('0')))
pclarv.ctd <- pclarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2==0 ~ as.numeric('0')))
pklarv.ctd <- pklarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2==0 ~ as.numeric('0')))
rxlarv.ctd <- rxlarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2==0 ~ as.numeric('0')))
yflarv.ctd <- yflarv.ctd %>% mutate(obs=case_when(Cper10m2>0 ~ as.numeric('1'),
                                                  Cper10m2==0 ~ as.numeric('0')))#making a new variable, "obs", that is binary for presence (1) or absence (0)


# Generalized Additive Model Formulations ---------------------------------

## Rex Sole, which had the worst performance of larval models:  
#models:
lv.base<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5),
             data=rxlarv.ctd,family=binomial)

lv.add.sal<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(salinity),
                data=rxlarv.ctd,family=binomial)

lv.add.temp<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+s(temperature),
                 data=rxlarv.ctd,family=binomial)

lv.2d<-gam(obs~factor(year)+s(doy)+s(lon,lat)+s(bottom_depth,k=5)+
             s(salinity,temperature),
           data=rxlarv.ctd,family=binomial)

#prediction: 
nlat=120
nlon=120
latd=seq(min(rxlarv.ctd$lat,na.rm=TRUE),max(rxlarv.ctd$lat,na.rm=TRUE),length.out=nlat)
lond=seq(min(rxlarv.ctd$lon,na.rm=TRUE),max(rxlarv.ctd$lon,na.rm=TRUE),length.out=nlon)

grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

grid.extent$dist<-NA
for(k in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[k],grid.extent$lon[k],
                          rxlarv.ctd$lat,rxlarv.ctd$lon)
  grid.extent$dist[k]<-min(dist)
}

grid.extent$year<-as.numeric(2016)
grid.extent$doy<-as.numeric(median(rxlarv.ctd$doy,na.rm=TRUE))
grid.extent$bottom_depth<-NA
grid.extent$bottom_depth<-as.numeric(median(rxlarv.ctd$bottom_depth,na.rm=TRUE))
grid.extent$pred<-predict(lv.base,newdata=grid.extent)
grid.extent$pred[grid.extent$dist>30000]<-NA

symcol<-adjustcolor('grey',alpha=0.5)

windows(height=15,width=15)
par(mai=c(1,1,0.5,0.9))
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),
                              ncol=length(lond),byrow=T)),col=hcl.colors(100,"PRGn"),
           ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),
           xlim=range(rxlarv.ctd$lon,na.rm=TRUE),ylim=range(rxlarv.ctd$lat,na.rm=TRUE),
           main='Alaska Plaice Larval Distribution, P/A Model',
           cex.main=1,cex.lab=1,cex.axis=0.9,legend.line=-2.5,
           legend.lab=expression(paste('Presence','- Absence')),legend.shrink=0.3)
contour(bathy,levels=-c(50,200),labcex=0.4,col='grey28',add=T)
points(rxlarv.ctd$lon[rxlarv.ctd$Cper10m2==0],rxlarv.ctd$lat[rxlarv.ctd$Cper10m2==0],pch='+',col='white')
symbols(rxlarv.ctd$lon[rxlarv.ctd$Cper10m2>0],
        rxlarv.ctd$lat[rxlarv.ctd$Cper10m2>0],
        circles=log(rxlarv.ctd$Cper10m2+1)[rxlarv.ctd$Cper10m2>0],
        inches=0.1,bg=symcol,fg='black',add=T)
map("worldHires",fill=T,col="seashell2",add=T)



















