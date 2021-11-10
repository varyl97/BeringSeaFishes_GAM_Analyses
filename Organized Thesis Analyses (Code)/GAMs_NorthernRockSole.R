# Generalized Additive Model Analyses for Northern Rock Sole:  -----------------------------------
#the following code creates generalized additive models for larvae of Pacific cod. 
#Pacific cod spawn beneath sea ice each year, therefore their eggs are rarely caught by ecoFOCI trawls. 
#   thus, Pac cod are only represented in the larval biogeography portion of my analyses. 

# Load larval data --------------------------------------------------------

nrslarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv',header=TRUE,check.names=TRUE)

##Now we create larval models, which are all additive but one has a two-dimensional smooth
#these models utilize CTD data
#begin with the larval base model: 
lv.base<-gam((Cper10m2+1)~factor(year)+s(doy)+s(lon,lat)+
               s(bottom_depth,k=5),
             data=nrslarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

windows(width=12,height=8)
par(mfrow=c(2,2))
plot(lv.base,select=1,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0,main='Base Larval Presence GAM')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.base,select=2,scheme=2,too.far=0.025,
     shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0)
map("world",fill=T,col="snow4",add=T)
plot(lv.base,select=3,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0)
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add in salinity in an additive interaction
lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                  s(bottom_depth,k=5)+
                  s(salinity),data=nrslarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

windows(width=12,height=8)
plot(lv.add.sal,page=1,scale=0,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE)

windows(width=12,height=8)
plot(lv.add.sal,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Salinity',
     xlab='Salinity (PSU)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add in temperature in an additive interaction
lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature),data=nrslarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

windows()
plot(lv.add.temp,page=1,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Larval Log(Cper10m2+1) w Temp',scale=0)
#note: scale=0 allows each individual smooth to have its own y-axis scale

windows(width=12,height=8)
plot(lv.add.temp,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Effect of Temperature',
     xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#Add both temperature and salinity in individual smooths each 
lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                   s(bottom_depth,k=5)+
                   s(temperature)+s(salinity),data=nrslarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

windows()
plot(lv.temp.sal,page=1,shade=TRUE,shade.col='skyblue4',
     main='Larval Log Presence Temp and Sal',
     seWithMean=TRUE,scale=0)

windows()
par(mfrow=c(1,2))
plot(lv.temp.sal,select=1,seWithMean=TRUE,shade=TRUE,
     shade.col='skyblue4',main='Temp and Sal Model')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
par(mfrow=c(1,2))
plot(lv.temp.sal,select=4,shade=TRUE,shade.col='skyblue4',
     seWithMean=TRUE,main='Effect of Temp',xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.temp.sal,select=5,shade=TRUE,shade.col='skyblue4',
     seWithMean=TRUE,main='Effect of Salinity',xlab='Salinity (psu)')
abline(h=0,col='sienna3',lty=2,lwd=2)

##Now try temperature and salinity in a two-dimensional, mutually interacting smooth: 
lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
             s(temperature,salinity),data=nrslarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)

windows()
par(mfrow=c(1,2))
plot(lv.2d,select=2,seWithMean=TRUE,shade=TRUE,shade.col='skyblue4',
     main='Seasonal Presence, 2D Temp+Sal Model')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.2d,select=1,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
plot(lv.2d,select=4,scheme=2,main='Larval Log Presence, 2D Temp and Sal Effect',
     too.far=0.025,
     xlab='Temperature (degC)',ylab='Salinity (psu)')


#checking based on AIC: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.tempsal<-AIC(lv.temp.sal)
aic.2d<-AIC(lv.2d)

aic.nrslarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                 'Add Sal, Temp','Sal-Temp 2D'),
                       'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                     aic.tempsal,aic.2d))

windows()
plot(c(1:5),aic.nrslarv$AIC_value,main='AIC Results for NRS Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.nrslarv$AIC_value,labels=round(aic.nrslarv$AIC_value),pos=c(4,1,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Add Sal, Temp','Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)


# Saved Larval Models -----------------------------------------------------
saveRDS(lv.base,"./GAM Models/nrs_larvae_base.rds")
saveRDS(lv.add.sal,"./GAM Models/nrs_larvae_addsal.rds")
saveRDS(lv.add.temp,"./GAM Models/nrs_larvae_addtemp.rds")
saveRDS(lv.temp.sal,"./GAM Models/nrs_larvae_addtempsal.rds")
saveRDS(lv.2d,"./GAM Models/nrs_larvae_2d.rds")

lv.base<-readRDS("./GAM Models/nrs_larvae_base.rds")
lv.add.sal<-readRDS("./GAM Models/nrs_larvae_addsal.rds")
lv.add.temp<-readRDS("./GAM Models/nrs_larvae_addtemp.rds")
lv.temp.sal<-readRDS("./GAM Models/nrs_larvae_addtempsal.rds")
lv.2d<-readRDS("./GAM Models/nrs_larvae_2d.rds")


