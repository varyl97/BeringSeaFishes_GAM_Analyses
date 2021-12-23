
# Kernel Density Analysis for Salinity-Temperature Plots ------------------

yflarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_YfLarv_wCTD.csv',header=TRUE,check.names=TRUE)
nrslarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_NrsLarv_wCTD.csv',header=TRUE,check.names=TRUE)
aplarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pclarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PcLarv_wCTD.csv',header=TRUE,check.names=TRUE)
pklarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_PkLarv_wCTD.csv',header=TRUE,check.names=TRUE)
fhlarv.ctd<-read.csv(file='./Ichthyo Data/Cleaned_Cut_FhLarv_wCTD.csv',header=TRUE,check.names=TRUE)

yfts<-ggplot(data=yflarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
  geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
  stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                  contour_var="ndensity",n=50,bins=15)+
  scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
  theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                   legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                   legend.justification = c(0, 1),
                   axis.ticks = element_line(colour = "grey70", size = 0.2),
                   panel.grid.major = element_line(colour = "grey70", size = 0.2),
                   panel.grid.minor = element_blank())+
  labs(x='Salinity',y='Temperature')+ggtitle("Yellowfin Sole")

nrsts<-ggplot(data=nrslarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
         geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
         stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                         contour_var="ndensity",n=50,bins=15)+
         scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
         theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                legend.justification = c(0, 1),
                axis.ticks = element_line(colour = "grey70", size = 0.2),
                panel.grid.major = element_line(colour = "grey70", size = 0.2),
                panel.grid.minor = element_blank())+
         labs(x='Salinity',y='Temperature')+ggtitle("Northern Rock Sole")

apts<-ggplot(data=aplarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
  geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
  stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                  contour_var="ndensity",n=50,bins=15)+
  scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
  theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                   legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                   legend.justification = c(0, 1),
                   axis.ticks = element_line(colour = "grey70", size = 0.2),
                   panel.grid.major = element_line(colour = "grey70", size = 0.2),
                   panel.grid.minor = element_blank())+
  labs(x='Salinity',y='Temperature')+ggtitle("Alaska Plaice")

pcts<-ggplot(data=pclarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
  geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
  stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                  contour_var="ndensity",n=50,bins=15)+
  scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
  theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                   legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                   legend.justification = c(0, 1),
                   axis.ticks = element_line(colour = "grey70", size = 0.2),
                   panel.grid.major = element_line(colour = "grey70", size = 0.2),
                   panel.grid.minor = element_blank())+
  labs(x='Salinity',y='Temperature')+ggtitle("Pacific Cod")

pkts<-ggplot(data=pklarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
  geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
  stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                  contour_var="ndensity",n=50,bins=15)+
  scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
  theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                   legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                   legend.justification = c(0, 1),
                   axis.ticks = element_line(colour = "grey70", size = 0.2),
                   panel.grid.major = element_line(colour = "grey70", size = 0.2),
                   panel.grid.minor = element_blank())+
  labs(x='Salinity',y='Temperature')+ggtitle("Walleye Pollock")

fhts<-ggplot(data=fhlarv.ctd,aes(x=salinity,y=temperature,weight=log(Cper10m2+1)))+
  geom_point(aes(size=log(Cper10m2+1)),color='darkgrey')+
  stat_density_2d(aes(fill=..level..),geom="polygon",color="white",alpha=0.3,
                  contour_var="ndensity",n=50,bins=15)+
  scale_fill_viridis_c(name='Density of log(Cper10m2+1)')+
  theme_bw()+theme(plot.title = element_text(face = "bold", size = 12),
                   legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                   legend.justification = c(0, 1),
                   axis.ticks = element_line(colour = "grey70", size = 0.2),
                   panel.grid.major = element_line(colour = "grey70", size = 0.2),
                   panel.grid.minor = element_blank())+
  labs(x='Salinity',y='Temperature')+ggtitle("Flathead Sole")


# Simplified 1-D Weighted KDE Estimates -----------------------------------

yfsal<-ggplot(data=yflarv.ctd)+geom_density(aes(x=salinity,weight=log(Cper10m2+1)))
yftemp<-ggplot(data=yflarv.ctd)+geom_density(aes(x=temperature,weight=log(Cper10m2+1)))
ggplotly(yfsal)
ggpltly(yftemp)


# Relic Code -------------------------------------------------------------------

test<-kde2d.weighted(nrslarv.ctd$salinity,nrslarv.ctd$temperature,w=log(nrslarv.ctd$Cper10m2+1))
test.df<-data.frame(expand.grid(x=test$x,y=test$y),z=as.vector(test$z))

windows()
ggplot(nrslarv.ctd,aes(x=salinity,y=temperature))+
  #geom_point()+
  geom_contour_filled(aes(x=x,y=y,z=z,fill=stat(level)),data=test.df)+
  xlim(29,35)+ylim(0,12)



# Box Plots for T,S Comparisons Across Species -----------------------------

yflarv.ctd$species<-as.character('YFS')
fhlarv.ctd$species<-as.character('FHS')
aplarv.ctd$species<-as.character('AP')
pklarv.ctd$species<-as.character('PK')
pclarv.ctd$species<-as.character('PC')
nrslarv.ctd$species<-as.character('NRS')

subyf<-yflarv.ctd[c('Cper10m2','salinity','temperature','species')]
subfh<-fhlarv.ctd[c('Cper10m2','salinity','temperature','species')]
subap<-aplarv.ctd[c('Cper10m2','salinity','temperature','species')]
subpk<-pklarv.ctd[c('Cper10m2','salinity','temperature','species')]
subpc<-pclarv.ctd[c('Cper10m2','salinity','temperature','species')]
subnrs<-nrslarv.ctd[c('Cper10m2','salinity','temperature','species')]

allspp<-rbind(subyf[subyf$Cper10m2>0,],subfh[subfh$Cper10m2>0,],
              subap[subap$Cper10m2>0,],subpk[subpk$Cper10m2>0,],
              subpc[subpc$Cper10m2>0,],subnrs[subnrs$Cper10m2>0,])

bpsal<-boxplot(allspp$salinity~allspp$species,xlab='Species',ylab='Salinity (psu)')
bpsal
windows()
par(oma=c(1,1.1,0.7,0.7))
boxplot(allspp$salinity~allspp$species,xlab='Species',ylab='Salinity (psu)')

bptemp<-boxplot(allspp$temperature~allspp$species,xlab='Species',
                ylab=expression(paste("Temperature ("^0, 'C)')))
bptemp
windows()
par(oma=c(1,1.1,0.7,0.7))
boxplot(allspp$temperature~allspp$species,xlab='Species',
        ylab=expression(paste("Temperature ("^0, 'C)')))


