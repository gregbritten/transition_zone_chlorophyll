rm(list=ls())

library(viridis)
library(maps)
library(fields)
library(lubridate)

source('r/load_data.r')
source('r/functions.r')

#########################################################################################
## AUGUST VS FEBRUARY ###################################################################
#########################################################################################
chllims <- c(-1.5,-0.0)
clims   <- c(-1,1.8)
cintlims<- c(-1,3.25)
chl_c_lims <- c(0,0.0275)
cols <- viridis(20)

##-GIOP-############################################
pdf('plots/seasonal_maps_giop.pdf',height=7,width=10,useDingbats=FALSE)
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
        lab2='d)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GIOP,2)),zlim=clims,
        lab1='',lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GIOP,2)/mapmean(CARBON_GIOP,2),zlim=chl_c_lims,
        lab1='',lab2='f)',lab3='')
mtext(outer=TRUE,side=1,'Longitude',line=0.5,cex=1.1)
mtext(outer=TRUE,side=2,'Latitude',line=0.5,cex=1.1)
dev.off()

##-GSM-#############################################
pdf('plots/seasonal_maps_gsm.pdf',height=7,width=10,useDingbats=FALSE)
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
        lab2='d)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CARBON_GSM,2)),zlim=clims,
        lab1='',lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=mapmean(CHL_GSM,2)/mapmean(CARBON_GSM,2),zlim=chl_c_lims,
        lab1='',lab2='f)',lab3='')
mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)
dev.off()


##--NPP MAPS--################################################
npplims <- c(1,3.5)
pdf('plots/seasonal_maps_npp.pdf',height=9,width=7,useDingbats=FALSE)
par(mfrow=c(4,2),mar=c(2,2,2,3.75),oma=c(3,3,2,2),xpd=FALSE)
plotmap(x=lon,y=lat,z=log10(mapmean(VGPM,8)),zlim=npplims,legend=FALSE,
        lab1=expression('log'['10']*'(VGPM NPP [mg/m'^2*'/day])'),
        lab2='a)',
        lab3='August')
plotmap(x=lon,y=lat,z=log10(mapmean(VGPM,2)),zlim=npplims,
        lab1=expression(''),
        lab2='b)',lab3='February')
plotmap(x=lon,y=lat,z=log10(mapmean(EP_VGPM,8)),zlim=npplims,legend=FALSE,
        lab1=expression('log'['10']*'(Eppley VGPM NPP [mg/m'^2*'/day])'),
        lab2='c)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(EP_VGPM,2)),zlim=npplims,legend=FALSE,
        lab1=expression(''),
        lab2='d)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CBPM,8)),zlim=npplims,legend=FALSE,
        lab1=expression('log'['10']*'(CBPM NPP [mg/m'^2*'/day])'),lab2='e)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CBPM,2)),zlim=npplims,legend=FALSE,
        lab1='',lab2='f)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CAFE,8)),zlim=npplims,legend=FALSE,
        lab1=expression('log'['10']*'(CAFE NPP [mg/m'^2*'/day])'),lab2='g)',lab3='')
plotmap(x=lon,y=lat,z=log10(mapmean(CAFE,2)),zlim=npplims,legend=FALSE,
        lab1='',lab2='h)',lab3='')


mtext(outer=TRUE,side=1,'Longitude',line=0.5)
mtext(outer=TRUE,side=2,'Latitude',line=0.5)
dev.off()


#######################################################################################
## CORRELATION MAPS ###################################################################
#######################################################################################
##--MAPS--###########################
pdf('plots/correlation_maps_giop_hycom.pdf',height=9,width=9)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plotcormap(cormap(CARBON_GIOP,MLD),lab='a) MLD vs. Carbon')
plotcormap(cormap(CHL_GIOP,MLD),lab='b) MLD vs. Chlorophyll')
plotcormap(cormap(CARBON_GIOP,I),lab='c) MLD-averaged Irradiance vs. Carbon')
plotcormap(cormap(CHL_GIOP,I),lab='c) MLD-averaged Irradiance vs. Chlorophyll')
mtext(side=1,outer=TRUE,line=0.5,'Longitude',cex=1.2)
mtext(side=2,outer=TRUE,line=0.5,'Latitude',cex=1.2)
dev.off()

pdf('plots/correlation_maps_giop_soda.pdf',height=9,width=9)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plotcormap(cormap(CARBON_GIOP[,,soda_i],ML_soda),lab='a) MLD vs. Carbon')
plotcormap(cormap(CHL_GIOP[,,soda_i],ML_soda),lab='b) MLD vs. Chlorophyll')
plotcormap(cormap(CARBON_GIOP[,,soda_i],I_soda),lab='c) MLD-averaged Irradiance vs. Carbon')
plotcormap(cormap(CHL_GIOP[,,soda_i],I_soda),lab='c) MLD-averaged Irradiance vs. Chlorophyll')
mtext(side=1,outer=TRUE,line=0.5,'Longitude',cex=1.2)
mtext(side=2,outer=TRUE,line=0.5,'Latitude',cex=1.2)
dev.off()

##--ZONAL AVERAGES--###############################
pdf('plots/zonal_correlations_giop_hycom.pdf',height=6,width=8)
par(mfrow=c(1,1),oma=c(2,2,2,2),cex.axis=0.75)
plot(lat,zonalcor(cormap(CHL_GIOP,MLD)),'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitude',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonalcor(cormap(CARBON_GIOP,MLD)), col='darkgreen')
lines(lat,zonalcor(cormap(CHL_GIOP,I)), col='green',lty=2)
lines(lat,zonalcor(cormap(CARBON_GIOP,I)), col='dark green', lty=2)
abline(h=0,lty=2)
legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. MLD-avg. Irradiance','Carbon vs. MLD-avg Irradiance'),
       col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2),cex=0.8)
dev.off()


pdf('plots/zonal_correlations_giop_soda.pdf',height=6,width=8)
par(mfrow=c(1,1),oma=c(2,2,2,2),cex.axis=0.75)
plot(lat,zonalcor(cormap(CHL_GIOP[,,soda_i],ML_soda)),'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitude',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonalcor(cormap(CARBON_GIOP[,,soda_i],ML_soda)), col='darkgreen')
lines(lat,zonalcor(cormap(CHL_GIOP[,,soda_i],I_soda)), col='green',lty=2)
lines(lat,zonalcor(cormap(CARBON_GIOP[,,soda_i],I_soda)), col='dark green', lty=2)
abline(h=0,lty=2)
legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. MLD-avg. Irradiance','Carbon vs. MLD-avg Irradiance'),
       col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2),cex=0.8)
dev.off()












