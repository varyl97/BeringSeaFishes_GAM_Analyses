########### Revised work with GAMs: 9/1/2021 ###########
world<-map_data("world")
BSmap<-world[world$long<(-155)&world$lat>50&world$lat<70,]
#syntax: +geom_map(data=BSmap,map=BSmap,aes(long,lat,map_id=region))

###EGGS: Spawning Behavior 
##Load in local and regional temperature index for March (2 mos before peak egg catch in May) 
loc.sst<-read.csv('../Environmental Data/Mar_SST_ByLocation_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(loc.sst) #more just to have, use regional for GAMs

reg.sst<-read.csv('../Environmental Data/Mar_SST_RegionalIndex_NCEP_BS.csv',header=TRUE,check.names=TRUE)
head(reg.sst) #range of regional average: lon: -180 to -151, lat: 50.5 to 67.5

for(i in 1:nrow(pksub)){
        pksub$reg.SST[i]<-reg.sst$SST[reg.sst$year==pksub$year[i]]}

pksub<-pksub[pksub$lat<60.5,]

#base: 
eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=pksub,family=tw(link='log'),method='REML')

summary(eg.base)

windows(width=12,height=8)
plot(eg.base,shade=TRUE,shade.col='skyblue3',page=1,
     seWithMean=TRUE,scale=0)


windows(width=12,height=8)#I think this is equivalent to plot above ^
vis.gam(eg.base,view=c(c("lon","lat")),too.far=0.025,
        color="heat",plot.type="contour",n.grid=20,
        main='Geographic Linear Predictor (vis.gam), Base Egg')
map("world",fill=T,col="snow4",add=T)

windows(width=12,height=8)
plot(eg.base,select=2,shade=TRUE,shade.col="skyblue3",
     xlab='Day of Year',ylab='Estimated Effect, Log(Cper10m2+1)',
     main="Phenologic Smooth, Base Egg Model")
abline(h=0,col='sienna3',lty=2,lwd=2)

windows()
par(mfrow=c(2,2),oma=c(1.5,2,1,1))
gam.check(eg.base)
mtext('gam.check Output for Base Egg PKE CPUE',outer=TRUE,cex=1,line=-0.5)

#different plots w/ Lorenzo's color gams code: 
temp<-tapply(pksub$reg.SST,list(pksub$year))

par(mfrow=c(5,5))

viridis<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF", "#25858EFF", "#1E9B8AFF",
          "#2BB07FFF", "#51C56AFF", "#85D54AFF", "#C2DF23FF", "#FDE725FF")) #viridis palette

windows(width=12,height=8)
myvis.gam(eg.base,view=c('lon','lat'),plot.type="contour",
          too.far=0.025,color='viridis',
          xlab='Longitude',
          ylab='Latitude',main='Base Egg Model')
map("world",fill=T,col="snow4",add=T)

plot(eg.base,select=2,xlab='Day of Year',shade=TRUE,shade.col='skyblue3')
abline(h=0,col='sienna3',lty=2,lwd=2)

#temp pheno threshold: 
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]
aic.pheno<-NA*(temps.in)
thr.pheno<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
        pksub$th<-factor(pksub$reg.SST<=temps.in[i])
        thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+
                                    s(lon,lat)+
                                    s(bottom_depth,k=5)+
                                    s(doy,by=th),
                                 data=pksub,family=tw(link='log'),method='REML')
        aic.pheno[i]<-AIC(thr.pheno[[i]])
}

best.index.phe<-order(aic.pheno)[1]

windows()
plot(temps.in,aic.pheno,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.pheno)),
     main='Temperature Threshold Flexible Phenology',ylab='AIC Index',
     xlab='Temperature (degC)')
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.phe],lty=2,lwd=2,col='steelblue3')

summary(thr.pheno[[best.index.phe]])

windows(width=12,height=8)
plot(thr.pheno[[best.index.phe]],shade=TRUE,shade.col='skyblue3',page=1,
     seWithMean=TRUE,scale=0)

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(thr.pheno[[best.index.phe]],select=4,main=paste('Below',temps.in[best.index.phe],sep=" "),
     shade=TRUE,shade.col='skyblue3',seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(thr.pheno[[best.index.phe]],select=3,main=paste('Above',temps.in[best.index.phe],sep=" "),
     shade=TRUE,shade.col='skyblue3',seWithMean=TRUE,xlab='Day of Year',ylab='Anomalies')
abline(h=0,col='sienna3',lty=2,lwd=2)

windows()
par(mfrow=c(2,2))
gam.check(thr.pheno[[best.index.phe]])

#temp threshold geo: 
temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]

aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
        pksub$th<-factor(pksub$reg.SST<=temps.in[i])
        thr.geo[[i]]<-gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                                  s(lon,lat,by=th),data=pksub,
                          family=tw(link='log'),method='REML')
        aic.geo[i]<-AIC(thr.geo[[i]])
}

best.index.geo<-order(aic.geo)[1]

windows()
plot(temps.in,aic.geo,type='b',lwd=2,ylim=range(c(AIC(eg.base),aic.geo)),
     main='Temperature Threshold Flex Geography',xlab="Temperature (degC)")
abline(h=AIC(eg.base),lty=2,lwd=2,col='sienna3')
abline(v=temps.in[best.index.geo],lty=2,lwd=2,col='steelblue3')

summary(thr.geo[[best.index.geo]])

windows(width=12,height=8)
plot(thr.geo[[best.index.geo]],page=1,scale=0,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE)

vis.gam(eg.base,view=c("lon","lat"),too.far=0.025,
        color="heat",plot.type="contour",n.grid=20,
        main='Geographic Linear Predictor (vis.gam), Base Egg')

windows(width=12,height=8)
par(mfrow=c(1,2))
plot(thr.geo[[best.index.geo]],select=4,scheme=2,too.far=0.025,
     main=paste('Below',temps.in[best.index.geo],sep=" "),
     shade=TRUE,seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)
plot(thr.geo[[best.index.geo]],select=3,scheme=2,too.far=0.025,
     main=paste('Above',temps.in[best.index.geo],sep=" "),
     shade=TRUE,seWithMean=TRUE,xlab='Longitude',ylab='Latitude')
map("world",fill=T,col="snow4",add=T)

windows()
par(mfrow=c(2,2))
gam.check(thr.geo[[best.index.geo]])

#vc temp pheno: 
vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                      s(doy,by=reg.SST),data=pksub,family=tw(link='log'),
              method='REML')
summary(vc.pheno)

windows(width=12,height=8)
plot(vc.pheno,shade=TRUE,shade.col='skyblue3',
     page=1,scale=0,main='V-C Temp Flexible Phenology, PK Eggs',
     seWithMean=TRUE)

windows(width=16,height=8)
par(mfrow=c(1,2))
plot(vc.pheno,select=2,shade=TRUE,shade.col='skyblue3',
     main='V-C Regional Temp, Flexible Phenology',xlab='Day of Year',
     seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(vc.pheno,select=4,shade=TRUE,shade.col='skyblue3',
     main='V-C Reg. Temp, Deviation from Avg. Pheno Variation',
     xlab='Day of Year',seWithMean=TRUE)
abline(h=0,col='sienna3',lty=2,lwd=2)

windows()
par(mfrow=c(2,2))
gam.check(vc.pheno)

#vc temp geo: 
vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                    s(lon,lat,by=reg.SST),data=pksub,family=tw(link='log'),
            method='REML')
summary(vc.geo)

windows(width=12,height=8)
plot(vc.geo,shade=TRUE,shade.col='skyblue3',
     page=1,seWithMean=TRUE,main='V-C Flex Geo')

windows(width=15,height=8)
par(mfrow=c(1,2))
plot(vc.geo,select=1,scheme=2,too.far=0.025,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,xlab='Longitude',ylab='Latitude',
     main='V-C Pk Egg Flex Geo, Avg. Variation')
map("world",fill=T,col="snow4",add=T)
plot(vc.geo,select=4,scheme=2,too.far=0.025,shade=TRUE,shade.col='skyblue3',
     xlab='Longitude',ylab='Latitude',seWithMean=TRUE,
     main='V-C Flex Geo, Deviation from Avg. Variation')
map("world",fill=T,col="snow4",add=T)

windows()
par(mfrow=c(2,2))
gam.check(vc.geo)

#all models for eggs in one place: 
eg.base<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5),
             data=pksub,family=tw(link='log'),method='REML')

summary(eg.base)

temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]
aic.pheno<-NA*(temps.in)
thr.pheno<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
        pksub$th<-factor(pksub$reg.SST<=temps.in[i])
        thr.pheno[[i]]<-gam((Cper10m2+1)~factor(year)+
                                    s(lon,lat)+
                                    s(bottom_depth,k=5)+
                                    s(doy,by=th),
                            data=pksub,family=tw(link='log'),method='REML')
        aic.pheno[i]<-AIC(thr.pheno[[i]])
}

best.index.phe<-order(aic.pheno)[1]
thr.pheno<-thr.pheno[[best.index.phe]]
summary(thr.pheno)

temps<-sort(unique(reg.sst$SST))
bd<-4
temps.in<-temps[bd:(length(temps)-bd)]

aic.geo<-NA*(temps.in)
thr.geo<-as.list(1:(length(temps.in)))

for(i in 1:length(temps.in)){
        pksub$th<-factor(pksub$reg.SST<=temps.in[i])
        thr.geo[[i]]<-gam((Cper10m2+1)~factor(year)+s(doy)+s(bottom_depth,k=5)+
                                  s(lon,lat,by=th),data=pksub,
                          family=tw(link='log'),method='REML')
        aic.geo[i]<-AIC(thr.geo[[i]])
}

best.index.geo<-order(aic.geo)[1]
thr.geo<-thr.geo[[best.index.geo]]
summary(thr.geo)

vc.pheno<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                      s(doy,by=reg.SST),data=pksub,family=tw(link='log'),
              method='REML')
summary(vc.pheno)

vc.geo<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy)+s(bottom_depth,k=5)+
                    s(lon,lat,by=reg.SST),data=pksub,family=tw(link='log'),
            method='REML')
summary(vc.geo)

#Save for reloading later!!
saveRDS(eg.base,file="../GAM Models/pk_egg_base.rds")
saveRDS(thr.pheno,file="../GAM Models/pk_egg_thr_pheno.rds")
saveRDS(thr.geo,file="../GAM Models/pk_egg_thr_geo.rds")
saveRDS(vc.pheno,file="../GAM Models/pk_egg_vc_pheno.rds")
saveRDS(vc.geo,file="../GAM Models/pk_egg_vc_geo.rds")

eg.base<-readRDS("../GAM Models/pk_egg_base.rds")
thr.pheno<-readRDS("../GAM Models/pk_egg_thr_pheno.rds")
thr.geo<-readRDS("../GAM Models/pk_egg_thr_geo.rds")
vc.pheno<-readRDS("../GAM Models/pk_egg_vc_pheno.rds")
vc.geo<-readRDS("../GAM Models/pk_egg_vc_geo.rds")

#checking based on AIC: 
aic.base<-AIC(eg.base)
aic.thrph<-AIC(thr.pheno[[best.index.phe]])
aic.thrge<-AIC(thr.geo[[best.index.geo]])
aic.vcph<-AIC(vc.pheno)
aic.vcgeo<-AIC(vc.geo)

aic.pkegg<-data.frame('model'=c('Base','Threshold Pheno','Threshold Geo',
                                'VC Pheno','VC Geo'),
                      'AIC_value'=c(aic.base,aic.thrph,aic.thrge,
                                    aic.vcph,aic.vcgeo))

windows()
plot(c(1:5),aic.pkegg$AIC_value,main='AIC Results for Pk Egg Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.pkegg$AIC_value,labels=round(aic.pkegg$AIC_value),pos=c(4,3,3,3,2))
legend("bottomleft",legend=c('Base','Threshold Pheno','Threshold Geo',
                              'VC Pheno','VC Geo'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)

###LARVAE: Water Mass Associations
lv.base<-gam((Cper10m2+1)~factor(year)+s(doy)+s(lon,lat)+
             s(bottom_depth,k=5),
             data=pklarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

windows(width=12,height=8)
par(mfrow=c(2,2))
plot(lv.base,select=1,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0,main='Base Larval Presence GAM, W/O Residuals')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.base,select=2,scheme=2,too.far=0.025,
     shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0)
map("world",fill=T,col="snow4",add=T)
plot(lv.base,select=3,shade=TRUE,shade.col='skyblue3',
     seWithMean=TRUE,scale=0)
abline(h=0,col='sienna3',lty=2,lwd=2)

#add salinity
lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                        s(bottom_depth,k=5)+
                        s(salinity),data=pklarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

windows(width=12,height=8)
plot(lv.add.sal,page=1,scale=0,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE)

windows(width=14,height=8)
par(mfrow=c(1,2))
plot(lv.add.sal,select=1,seWithMean=TRUE,shade=TRUE,shade.col="skyblue4",
     main='Seasonal Presence, Added Sal Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.add.sal,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows(width=12,height=8)
plot(lv.add.sal,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Salinity',
     xlab='Salinity (PSU)')
abline(h=0,col='sienna3',lty=2,lwd=2)

windows(width=12,height=8)
plot(lv.add.sal,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Effect of Salinity, W/O Residuals',
     xlab='Salinity (PSU)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#add temperature
lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                         s(bottom_depth,k=5)+
                         s(temperature),data=pklarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

windows()
plot(lv.add.temp,page=1,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Larval Log(Cper10m2+1) w Temp',scale=0)

windows()
plot(lv.add.temp,page=1,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='W/O Residuals',scale=0)

windows(width=14,height=8)
par(mfrow=c(1,2))
plot(lv.add.temp,select=1,seWithMean=TRUE,shade=TRUE,shade.col="skyblue4",
     main='Seasonal Presence, Added Temp Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.add.temp,select=2,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows(width=12,height=8)
plot(lv.add.temp,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Larval Log(CPer10m2+1), Effect of Temperature',
     xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)

windows(width=12,height=8)
plot(lv.add.temp,select=4,shade=TRUE,shade.col="skyblue4",
     seWithMean=TRUE,main='Effect of Temperature, W/O Residuals',
     xlab='Temperature (degC)')
abline(h=0,col='sienna3',lty=2,lwd=2)

#additive both temp and sal
lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                         s(bottom_depth,k=5)+
                         s(temperature)+s(salinity),data=pklarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

windows()
plot(lv.temp.sal,page=1,shade=TRUE,shade.col='skyblue4',
     main='Larval Log Presence Temp and Sal',
     seWithMean=TRUE,scale=0)

windows()
plot(lv.temp.sal,page=1,shade=TRUE,shade.col='skyblue4',
     main='Temp and Sal, w/o Resids',seWithMean=TRUE,scale=0)

windows()
par(mfrow=c(1,2))
plot(lv.temp.sal,select=1,seWithMean=TRUE,shade=TRUE,
     shade.col='skyblue4',main='Temp and Sal Model (No resids.)')
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

##2D Smooth with temp and sal: 
lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
                   s(temperature,salinity),data=pklarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)

windows()
plot(lv.2d,page=1,shade=TRUE,shade.col='skyblue4',
     main='Larval Log Presence, 2D Temp and Sal',
     seWithMean=TRUE,scale=0)

windows()
par(mfrow=c(1,2))
plot(lv.2d,select=2,seWithMean=TRUE,shade=TRUE,shade.col='skyblue4',
     main='Seasonal Presence, 2D Temp+Sal Model (No resids.)')
abline(h=0,col='sienna3',lty=2,lwd=2)
plot(lv.2d,select=1,scheme=2,seWithMean=TRUE,too.far=0.025,
     xlab='Longitude',ylab='Latitude',main='Biogeography')
map("world",fill=T,col="snow4",add=T)

windows()
plot(lv.2d,select=4,scheme=2,main='Larval Log Presence, 2D Temp and Sal Effect',
     too.far=0.025,
     xlab='Temperature (degC)',ylab='Salinity (psu)')

#all larval models in one place: 
lv.base<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                     s(bottom_depth,k=5),
             data=pklarv.ctd,family=tw(link='log'),method='REML')
summary(lv.base)

lv.add.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                        s(bottom_depth,k=5)+
                        s(salinity),data=pklarv.ctd,family=tw(link='log'),
                method='REML')
summary(lv.add.sal)

lv.add.temp<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                         s(bottom_depth,k=5)+
                         s(temperature),data=pklarv.ctd,family=tw(link='log'),
                 method='REML')
summary(lv.add.temp)

lv.temp.sal<-gam((Cper10m2+1)~factor(year)+s(doy,k=7)+s(lon,lat)+
                         s(bottom_depth,k=5)+
                         s(temperature)+s(salinity),data=pklarv.ctd,
                 family=tw(link='log'),method='REML')
summary(lv.temp.sal)

lv.2d<-gam((Cper10m2+1)~factor(year)+s(lon,lat)+s(doy,k=7)+s(bottom_depth)+
                   s(temperature,salinity),data=pklarv.ctd,family=tw(link='log'),
           method='REML')
summary(lv.2d)

#checking based on AIC: 
aic.base.lv<-AIC(lv.base)
aic.sal<-AIC(lv.add.sal)
aic.temp<-AIC(lv.add.temp)
aic.tempsal<-AIC(lv.temp.sal)
aic.2d<-AIC(lv.2d)

aic.pklarv<-data.frame('model'=c('Base BioGeo','Add Sal','Add Temp',
                                'Add Sal, Temp','Sal-Temp 2D'),
                      'AIC_value'=c(aic.base.lv,aic.sal,aic.temp,
                                    aic.tempsal,aic.2d))

windows()
plot(c(1:5),aic.pklarv$AIC_value,main='AIC Results for Pk Larvae Models',
     col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
     pch=19,cex=2,ylab='AIC Value',xlab='')
grid(nx=5,ny=14,col="lightgray")
text(c(1:5),aic.pklarv$AIC_value,labels=round(aic.pklarv$AIC_value),pos=c(4,3,3,3,2))
legend("bottomleft",legend=c('Base BioGeo','Add Sal','Add Temp',
                             'Add Sal, Temp','Sal-Temp 2D'),
       col=c( "#482173FF", "#38598CFF","#1E9B8AFF", "#51C56AFF","#FDE725FF"),
       lwd=3,lty=1)

#finding salinity/temperature hotspots 
viridis<-colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF", "#2D708EFF", "#25858EFF", "#1E9B8AFF",
                            "#2BB07FFF", "#51C56AFF", "#85D54AFF", "#C2DF23FF", "#FDE725FF")) #viridis palette

stsweet<-pklarv.ctd[pklarv.ctd$temperature>4.4&pklarv.ctd$temperature<5.6&
                            pklarv.ctd$salinity>31.4&pklarv.ctd$salinity<32.6,]#dim 249, 27
tt<-table(stsweet$year)

stsweet<-subset(stsweet,year%in%names(tt[tt>5]))
table(stsweet$year) #only includes 2002,2003,2005,2009,2012,2014,2015,2016 
stsweet<-subset(stsweet,Cper10m2>0)

stsweet$color<-NA
stsweet$month_nm<-NA
stsweet<-stsweet %>% mutate(color=case_when(
        stsweet$month==3~'#440154FF',
        stsweet$month==4~'#2D708EFF',
        stsweet$month==5~'#1E9B8AFF',
        stsweet$month==6~'#C2DF23FF'))
stsweet<-stsweet %>%mutate(month_nm=case_when(
        stsweet$month==3~'March',
        stsweet$month==4~'April',
        stsweet$month==5~'May',
        stsweet$month==6~'June'))

test<-stsweet[stsweet$year==2002,]

windows()
plot(test$lon,test$lat,xlim=range(pklarv.ctd$lon),ylim=range(pklarv.ctd$lat),
     col=test$color,pch=19,cex=2,
     xlab='Longitude',ylab='Latitude',main='Test, Salinity/Temperature Hot Spots')
symbols(test$lon,test$lat,circles=test$Cper10m2,inches=0.15,add=T)
map("worldHires",fill=T,col="snow4",add=T) 
#NOTE: when bathymetry is figured out, can add station points via points()

years<-sort(unique(stsweet$year))
tmp1<-1:ceiling(length(years)/4)
lg.text<-unique(stsweet$month_nm)
lg.values<-c('#440154FF','#2D708EFF','#1E9B8AFF',
             '#C2DF23FF')

for(j in 1:length(tmp1)){
        windows(width=24,height=14)
        par(mfcol=c(2,2),omi=c(0.25,0.3,0.55,0.25),mai=c(0.2,0.4,0.4,0.1))
        for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
         plot(stsweet$lon[stsweet$year==years[i]],stsweet$lat[stsweet$year==years[i]],
              col=stsweet$color[stsweet$year==years[i]],pch=19,cex=2,
              main=as.character(years[i]),
              ylim=range(stsweet$lat),xlim=range(stsweet$lon),
              ylab=expression(paste("Latitude ("^0,'N)')),
              xlab=expression(paste("Longitude ("^0,'E)')))
         symbols(pklarv.ctd$lon[pklarv.ctd$Cper10m2>0&pklarv.ctd$year==years[i]],
                 pklarv.ctd$lat[pklarv.ctd$Cper10m2>0&pklarv.ctd$year==years[i]],
                 circles=pklarv.ctd$Cper10m2[pklarv.ctd$Cper10m2>0&pklarv.ctd$year==years[i]],
                 inches=0.15,add=T)
         map("worldHires",fill=T,col="snow4",add=T)
        }
    legend("bottomright",legend=c('March','April','May','June'),
               col=lg.values,lwd=3,lty=1,
               bg='lightgrey')
    mtext("Pollock +Cper10m2 Catch Plotted Against T-S 'Sweet Spots'",outer=TRUE,cex=1,line=1)
} 

























