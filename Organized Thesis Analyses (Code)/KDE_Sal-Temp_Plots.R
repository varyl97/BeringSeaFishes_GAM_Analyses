
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


# Relic Code -------------------------------------------------------------------

test<-kde2d.weighted(nrslarv.ctd$salinity,nrslarv.ctd$temperature,w=log(nrslarv.ctd$Cper10m2+1))
test.df<-data.frame(expand.grid(x=test$x,y=test$y),z=as.vector(test$z))

windows()
ggplot(nrslarv.ctd,aes(x=salinity,y=temperature))+
  #geom_point()+
  geom_contour_filled(aes(x=x,y=y,z=z,fill=stat(level)),data=test.df)+
  xlim(29,35)+ylim(0,12)



