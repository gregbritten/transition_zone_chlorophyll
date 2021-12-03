rm(list=ls())

library(Cairo)
library(lubridate)
library(viridis)
library(maps)
library(fields)

source('r/load_data.r')
source('r/functions.r')

chl_feb  <- CHL[,,which(month(date)==2)]
chl_aug  <- CHL[,,which(month(date)==8)]
c_feb    <- CARBON[,,which(month(date)==2)]
c_aug    <- CARBON[,,which(month(date)==8)]

lat <- seq(0,65,length.out=130)
lon <- seq(-180,-115,length.out=130)

chl_aug_mean <- apply(chl_aug,c(1,2),function(x) mean(x,na.rm=TRUE))
chl_feb_mean <- apply(chl_feb,c(1,2),function(x) mean(x,na.rm=TRUE))
c_aug_mean <- apply(c_aug,c(1,2),function(x) mean(x,na.rm=TRUE))
c_feb_mean <- apply(c_feb,c(1,2),function(x) mean(x,na.rm=TRUE))

chl_c_aug <- chl_aug_mean/c_aug_mean
chl_c_feb <- chl_feb_mean/c_feb_mean

#########################################################################################
## FIGURE 1 #############################################################################
#########################################################################################
chllims <- c(-1.5,-0.0)
clims   <- c(-2,1.6)
cintlims<- c(-1,3.25)
chl_c_lims <- c(0,0.0275)
cols <- viridis(20)

pdf('~/dropbox/working/gradients/tzcf/plots/maps_05_19_2021.pdf',height=7,width=10,useDingbats=FALSE)
par(mfrow=c(2,3),mar=c(2,2,2,5),oma=c(3,3,2,2),xpd=FALSE)
image2(x=lon,y=lat,z=log10(chl_aug_mean),col=viridis(25),zlim=chllims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  box()
  mtext(adj=0,expression('log'['10']*'(Chlorophyll [mg/m'^3*'])'),cex=0.8)
  mtext(adj=-0.1,'a)')
  mtext(adj=0,'August',line=1.5)
  image.plot(matrix(chllims),legend.only=TRUE,col=cols)
image2(x=lon,y=lat,z=log10(c_aug_mean),col=viridis(25),zlim=clims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  box()
  mtext(adj=0,expression('log'['10']*'(Carbon Concentration [mg/m'^3*'])'),cex=0.8)
  mtext(adj=-0.1,'b)')
  image.plot(matrix(clims),legend.only=TRUE,col=cols)
image2(x=lon,y=lat,z=chl_c_aug,col=viridis(25),zlim=chl_c_lims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  box()
  mtext(adj=0,expression('Chlorophyll:Carbon'),cex=0.8)
  mtext(adj=-0.1,'c)')
  image.plot(matrix(chl_c_lims),legend.only=TRUE,col=cols)

image2(x=lon,y=lat,z=log10(chl_feb_mean),col=viridis(25),zlim=chllims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  box()
  mtext(adj=-0.15,expression('d)'))
  mtext(adj=0,'February',line=0.4)
image2(x=lon,y=lat,z=log10(c_feb_mean),col=viridis(25),zlim=clims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  mtext(adj=-0.1,expression('e)'))Áˇ
  box()
image2(x=lon,y=lat,z=chl_c_feb,col=viridis(25),zlim=chl_c_lims)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  mtext(adj=-0.1,expression('f)'))
  box()

mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)

dev.off()


