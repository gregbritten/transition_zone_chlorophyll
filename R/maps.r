rm(list=ls())

library(Cairo)
library(lubridate)
library(viridis)
library(maps)
library(fields)

image2 <- function(x,y,z,zlim,col){
  tmp <- z
  tmp[tmp>zlim[2]] <- zlim[2]
  tmp[tmp<zlim[1]] <- zlim[1]
  image(x,y,tmp,col=col)
}


#source('d:/dropbox/working/gradients/tzcf/github/greg/load_data.r')
source('~/dropbox/working/gradients/tzcf/github/greg/load_data.r')
load('~/dropbox/working/gradients/tzcf/data/CBPM_global.rdata')
load('~/dropbox/working/gradients/tzcf/data/VGPM_global.rdata')


########################################################################
## COMPARE CBPM AND VGPM ###############################################
########################################################################
CBPM_global_clim <- apply(CBPM_global,c(1,2),function(x) mean(x,na.rm=TRUE))
VGPM_global_clim <- apply(VGPM_global,c(1,2),function(x) mean(x,na.rm=TRUE))

pdf('~/dropbox/working/gradients/tzcf/plots/cbpm_vgpm_climatology.pdf',height=9,width=10)
par(mfrow=c(2,1),mar=c(2,2,2,4),oma=c(2,2,2,8))
image2(y=seq(-90,90,length.out=360),x=seq(-180,180,length.out=360*2),
       z=apply(log10(CBPM_global_clim),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(col='grey',fill=TRUE,add=TRUE) 
box()
mtext(adj=0,expression('log'['10']*'(CBPM NPP [mgC/m'^2*'/day])'))

image.plot(matrix(npplims),legend.only=TRUE,col=cols)
image2(y=seq(-90,90,length.out=360),x=seq(-180,180,length.out=360*2),
       z=apply(log10(VGPM_global_clim),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(col='grey',fill=TRUE,add=TRUE) 
box()
mtext(adj=0,expression('log'['10']*'(VGPM NPP [mgC/m'^2*'/day])'))

dev.off()



chl_feb  <- CHL[,,which(month(ocnprd_date)==2)]
chl_aug  <- CHL[,,which(month(ocnprd_date)==8)]
c_feb    <- C[,,which(month(ocnprd_date)==2)]
c_aug    <- C[,,which(month(ocnprd_date)==8)]
cint_feb <- Cint[,,which(month(ocnprd_date)==2)]
cint_aug <- Cint[,,which(month(ocnprd_date)==8)]
cbpm_feb <- CBPM[,,which(month(ocnprd_date)==2)]
cbpm_aug <- CBPM[,,which(month(ocnprd_date)==8)]
cafe_feb <- CAFE[,,which(month(ocnprd_date)==2)]
cafe_aug <- CAFE[,,which(month(ocnprd_date)==8)]
vgpm_feb <- VGPM[,,which(month(ocnprd_date)==2)]
vgpm_aug <- VGPM[,,which(month(ocnprd_date)==8)]
evgpm_feb <- EPPLEY_VGPM[,,which(month(ocnprd_date)==2)]
evgpm_aug <- EPPLEY_VGPM[,,which(month(ocnprd_date)==8)]
bbp_feb <- BBP[,,which(month(ocnprd_date)==2)]
bbp_aug <- BBP[,,which(month(ocnprd_date)==8)]

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
#npplims <- c(1.0,3.5)
#bbplims <- c(0,0.025)
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


sudo ln -s /Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia


#######################################################
## NPP MAPS ###########################################
#######################################################
pdf('~/dropbox/working/gradients/tzcf/plots/maps_npp_03_26_2021.pdf',height=10,width=7,useDingbats=FALSE)
par(mfrow=c(4,2),mar=c(2,2,2,4),oma=c(2,2,2,2))

image2(x=lon,y=lat,z=apply(log10(vgpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
mtext(adj=0,expression('log'['10']*'(VGPM NPP [mgC/m'^2*'/day])'))
mtext(adj=0,'August',line=1.5)
image2(x=lon,y=lat,z=apply(log10(vgpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
image.plot(matrix(npplims),legend.only=TRUE,col=cols)
mtext(adj=0,'February',line=1.5)

image2(x=lon,y=lat,z=apply(log10(evgpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
mtext(adj=0,expression('log'['10']*'(Eppley VGPM NPP [mgC/m'^2*'/day])'))
image2(x=lon,y=lat,z=apply(log10(evgpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
#image.plot(matrix(npplims),legend.only=TRUE,col=cols)


image2(x=lon,y=lat,z=apply(log10(cbpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
mtext(adj=0,expression('log'['10']*'(CBPM NPP [mgC/m'^2*'/day])'))
image2(x=lon,y=lat,z=apply(log10(cbpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
#image.plot(matrix(npplims),legend.only=TRUE,col=cols)
#
image2(x=lon,y=lat,z=apply(log10(cafe_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
mtext(adj=0,expression('log'['10']*'(CAFE NPP [mgC/m'^2*'/day])'))
image2(x=lon,y=lat,z=apply(log10(cafe_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE)
box()
#image.plot(matrix(npplims),legend.only=TRUE,col=cols)

dev.off()










####################################################
#CairoPDF('d:/dropbox/working/gradients/tzcf/plots/maps_npp_01_28_2021.pdf',height=7,width=6)
# pdf('~/dropbox/working/gradients/tzcf/plots/maps_npp_01_28_2021.pdf',height=7,width=6,useDingbats=FALSE)
# par(mfrow=c(3,2),mar=c(2,2,2,4),oma=c(2,2,2,2))
# image2(x=lon,y=lat,z=apply(log10(cbpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# mtext(adj=0,expression('log'['10']*'(CBPM NPP [mgC/m'^2*'/day])'))
# image2(x=lon,y=lat,z=apply(log10(cbpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# image.plot(matrix(npplims),legend.only=TRUE,col=cols)

# image2(x=lon,y=lat,z=apply(log10(cafe_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# mtext(adj=0,expression('log'['10']*'(CAFE NPP [mgC/m'^2*'/day])'))
# image2(x=lon,y=lat,z=apply(log10(cafe_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# image.plot(matrix(npplims),legend.only=TRUE,col=cols)

# image2(x=lon,y=lat,z=apply(log10(vgpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# mtext(adj=0,expression('log'['10']*'(VGPM NPP [mgC/m'^2*'/day])'))
# image2(x=lon,y=lat,z=apply(log10(vgpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
# map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
# box()
# image.plot(matrix(npplims),legend.only=TRUE,col=cols)

#image2(x=lon,y=lat,z=apply(log10(evgpm_aug),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
#map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
#box()
#mtext(adj=0,expression('log'['10']*'(VGPM NPP [mgC/m'^2*'/day])'))
#image2(x=lon,y=lat,z=apply(log10(evgpm_feb),c(1,2),function(x) mean(x,na.rm=TRUE)),col=viridis(25),zlim=npplims)
#map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
#box()
#image.plot(matrix(npplims),legend.only=TRUE,col=cols)

# mtext(outer=TRUE,adj=0.2,'August')
# mtext(outer=TRUE,adj=0.75,'February')
# 
# dev.off()

