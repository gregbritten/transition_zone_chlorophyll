rm(list=ls())

library(Cairo)
library(lubridate)
library(viridis)
library(maps)
library(fields)

source('r/load_data.r')
source('r/functions.r')

#########################################################################################
## FIGURE 1 #############################################################################
#########################################################################################
chllims <- c(-1.5,-0.0)
clims   <- c(-0,2)
cintlims<- c(-1,3.25)
chl_c_lims <- c(0,0.0275)
cols <- viridis(20)

pdf('~plots/maps.pdf',height=7,width=10,useDingbats=FALSE)
par(mfrow=c(2,3),mar=c(2,2,2,5),oma=c(3,3,2,2),xpd=FALSE)
plotmap(x=lon,y=lat,z=log10(mapmean(CHL,8)),zlim=chllims,
        lab1=expression('log'['10']*'(Chlorophyll [mg/m'^3*'])'),
        lab2='a)',
        lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON,8)),zlim=clims,
        lab1=expression('log'['10']*'(Carbon Concentration [mg/m'^3*'])'),
        lab2='b)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL,8)/mapmean(CARBON,8),zlim=chl_c_lims,
        lab1=expression('Chlorophyll:Carbon'),
        lab2='c)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CHL,2)),zlim=chllims,
        lab1=expression('February'),
        lab2='d)',lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON,2)),zlim=clims,
        lab1='',lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL,2)/mapmean(CARBON,2),zlim=chl_c_lims,
        lab1='',lab2='f)',lab3='')
mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)
dev.off()

####################################################
## GIOP ############################################
####################################################
par(mfrow=c(2,3),mar=c(2,2,2,5),oma=c(3,3,2,2),xpd=FALSE)
plotmap(x=lon,y=lat,z=log10(mapmean(CHL_GIOP,8)),zlim=chllims,
        lab1=expression('log'['10']*'(Chlorophyll [mg/m'^3*'])'),
        lab2='a)',
        lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GIOP,8)),zlim=clims,
        lab1=expression('log'['10']*'(Carbon Concentration [mg/m'^3*'])'),
        lab2='b)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GIOP,8)/mapmean(CARBON_GIOP,8),zlim=chl_c_lims,
        lab1=expression('Chlorophyll:Carbon'),
        lab2='c)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CHL_GIOP,2)),zlim=chllims,
        lab1=expression('February'),
        lab2='d)',lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GIOP,2)),zlim=clims,
        lab1='',lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GIOP,2)/mapmean(CARBON_GIOP,2),zlim=chl_c_lims,
        lab1='',lab2='f)',lab3='')
mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)

####################################################
## GSM #############################################
####################################################
par(mfrow=c(2,3),mar=c(2,2,2,5),oma=c(3,3,2,2),xpd=FALSE)
plotmap(x=lon,y=lat,z=log10(mapmean(CHL_GSM,8)),zlim=chllims,
        lab1=expression('log'['10']*'(Chlorophyll [mg/m'^3*'])'),
        lab2='a)',
        lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GSM,8)),zlim=clims,
        lab1=expression('log'['10']*'(Carbon Concentration [mg/m'^3*'])'),
        lab2='b)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GSM,8)/mapmean(CARBON_GSM,8),zlim=chl_c_lims,
        lab1=expression('Chlorophyll:Carbon'),
        lab2='c)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CHL_GSM,2)),zlim=chllims,
        lab1=expression('February'),
        lab2='d)',lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GSM,2)),zlim=clims,
        lab1='',lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GSM,2)/mapmean(CARBON_GSM,2),zlim=chl_c_lims,
        lab1='',lab2='f)',lab3='')
mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)














