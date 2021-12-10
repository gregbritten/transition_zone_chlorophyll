rm(list=ls())
library(viridis)
library(fields)
library(RColorBrewer)

source("r/load_data.r")
source("r/functions.r")

labs <- c('f) 00-10 Deg. N', 'e) 10-20 Deg. N', 'd) 20-30 Deg. N', 'c) 30-40 Deg. N', 'b) 40-50 Deg. N', 'a) 50-60 Deg. N')
timedims <- c(2003,2021)
ticks <- 2003:2021

######################################################################
## COMPARE GIOP and GSM ##############################################
######################################################################
##--CHLOROPHYLL
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  chl_gsmi   <- latmean(CHL_GSM, loni=loni,lati=LAT[[i]])
  chl_giopi  <- latmean(CHL_GIOP,loni=loni,lati=LAT[[i]])
  plot(time,chl_gsmi,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n',ylim=range(c(chl_gsmi,chl_giopi),na.rm=TRUE)); 
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=4)
  lines(time,chl_giopi,type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chl_gsmi,chl_giopi,use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
  axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))

##--CARBON
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  c_gsmi   <- latmean(CARBON_GSM, loni=loni,lati=LAT[[i]])
  c_giopi  <- latmean(CARBON_GIOP,loni=loni,lati=LAT[[i]])
  plot(time,c_gsmi,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n',ylim=range(c(c_gsmi,c_giopi),na.rm=TRUE)); 
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=4)
  lines(time,c_giopi,type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(c_gsmi,c_giopi,use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
  axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))

##--RATIO
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  chl_gsmi   <- latmean(CHL_GSM, loni=loni,lati=LAT[[i]])
  chl_giopi  <- latmean(CHL_GIOP,loni=loni,lati=LAT[[i]])
  c_gsmi     <- latmean(CARBON_GSM, loni=loni,lati=LAT[[i]])
  c_giopi    <- latmean(CARBON_GIOP,loni=loni,lati=LAT[[i]])
  plot(time,chl_gsmi/c_gsmi,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n',ylim=range(c(chl_gsmi/c_gsmi,chl_giopi/c_giopi),na.rm=TRUE)); 
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=4)
  lines(time,chl_giopi/c_giopi,type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chl_gsmi/c_gsmi,chl_giopi/c_giopi,use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
  axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))

#####################################################################
## COMPARISON OF NPP PRODUCTS #######################################
#####################################################################
cols <- brewer.pal(4,"Dark2")
pdf('plots/chl_npp_comparison.pdf',height=10,width=8)
par(mfrow=c(12,1),mar=c(1,2,1,2),oma=c(2,3,2,5),cex.axis=0.7)
for(i in 6:1){
  cbpmi   <- latmean(CBPM,loni=loni,lati=LAT[[i]])
  vgpmi   <- latmean(VGPM,loni=loni,lati=LAT[[i]])
  cafei   <- latmean(CAFE,loni=loni,lati=LAT[[i]])
  epvgpmi <- latmean(VGPM,loni=loni,lati=LAT[[i]])
  chli    <- latmean(CHL,loni=loni,lati=LAT[[i]])
  plot(time,cbpmi,type='l',col=cols[1],xlim=timedims,xaxt='n',yaxt='n',ylim=range(c(cbpmi,vgpmi,cafei,epvgpmi),na.rm=TRUE)); 
  lines(time,vgpmi,type='l',col=cols[2],xaxt='n',xlim=timedims,yaxt='n')
  lines(time,cafei,type='l',col=cols[3],xaxt='n',xlim=timedims,yaxt='n')
  lines(time,epvgpmi,type='l',col=cols[4],xaxt='n',xlim=timedims,yaxt='n')
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=2)
  par(new=TRUE)	
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chli,cbpmi,use='pairwise.complete.obs'),3)),adj=0.7,cex=0.8,col=cols[1])	
  mtext(paste(round(cor(chli,vgpmi,use='pairwise.complete.obs'),3)),adj=0.8,cex=0.8,col=cols[2])	
  mtext(paste(round(cor(chli,cafei,use='pairwise.complete.obs'),3)),adj=0.9,cex=0.8,col=cols[3])	
  mtext(paste(round(cor(chli,epvgpmi,use='pairwise.complete.obs'),3)),adj=1,cex=0.8,col=cols[4])	
plot(time,chli,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
axis(1,labels=NA,at=ticks)
abline(v=ticks,lty=2)
axis(side=4)
}
axis(1,at=ticks)
mtext(side=4,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'),line=0.75)
mtext(side=2,outer=TRUE, expression('Net Primary Productivity [mg/m'^2*'/day]'))
dev.off()


#####################################################################
## TIME SERIES PLOTS - CHL VS NPP CBPM ##############################
#####################################################################
pdf('plots/chl_npp.pdf',height=6,width=5)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  cbpmi <- latmean(CBPM,loni=loni,lati=LAT[[i]])
  chli  <- latmean(CHL,loni=loni,lati=LAT[[i]])
  plot(time,cbpmi,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
	abline(v=ticks,lty=2)
	axis(1,labels=NA,at=ticks)
	axis(side=4)
par(new=TRUE)	
plot(time,chli,type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
	mtext(adj=0,labs[i],cex=0.6)
	axis(1,labels=NA,at=ticks)
	mtext(paste(round(cor(chli,cbpmi,use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
	axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(side=4,outer=TRUE, expression('Net Primary Productivity [mg/m'^2*'/day]'),line=0.75)
dev.off()

#####################################################################
## TIME SERIES PLOTS - CHL VS NPP VGPM ##############################
#####################################################################
pdf('~/dropbox/working/gradients/tzcf/plots/chl_npp_vgpm.pdf',height=6,width=5)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  vgpmi <- latmean(VGPM,loni=loni,lati=LAT[[i]])
  chli  <- latmean(CHL,loni=loni,lati=LAT[[i]])
  plot(time,vgpmi,type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=4)
  par(new=TRUE)	
  plot(time,chli,type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chli,vgpmi,use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
  axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(side=4,outer=TRUE, expression('Net Primary Productivity [mg/m'^2*'/day]'),line=0.75)
dev.off()

#####################################################################
## TIME SERIES PLOTS - CHL VS CARBON ################################
#####################################################################
chllims  <- c(0,0.4)
timedims <- c(2003,2020)
cdims    <- c(5,30) 

pdf('~plots/c_vs_chl_timeseries.pdf',height=6,width=8)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  plot(time,latmean(CHL,loni=loni,lati=LAT[[i]]),type='l',col='dark green',yaxt='n',xlim=timedims,xaxt='n'); 
  axis(side=2,col='dark green'); 
  axis(1,labels=NA,at=ticks)
  
  par(new=TRUE)
  plot(time,latmean(CARBON,loni=loni,lati=LAT[[i]]),type='l',col='black',xaxt='n',yaxt='n',xlim=timedims)
  axis(4,col='black')
  abline(v=ticks,lty=2)
  mtext(adj=0,labs[i],cex=0.6)
  mtext(paste(round(cor(chl[[i]],carbon[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)
}
axis(1,at=ticks)
mtext(outer=TRUE, side=2, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(outer=TRUE, side=4, expression('Carbon Concentration [mg/m'^3*']'),line=0.75)
dev.off()

##--GIOP--################
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  plot(time,latmean(CHL_GIOP,loni=loni,lati=LAT[[i]]),type='l',col='dark green',yaxt='n',xlim=timedims,xaxt='n'); 
  axis(side=2,col='dark green'); 
  axis(1,labels=NA,at=ticks)
  
  par(new=TRUE)
  plot(time,latmean(CARBON_GIOP,loni=loni,lati=LAT[[i]]),type='l',col='black',xaxt='n',yaxt='n',xlim=timedims)
  axis(4,col='black')
  abline(v=ticks,lty=2)
  mtext(adj=0,labs[i],cex=0.6)
  mtext(paste(round(cor(chl[[i]],carbon[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)
}
axis(1,at=ticks)
mtext(outer=TRUE, side=2, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(outer=TRUE, side=4, expression('Carbon Concentration [mg/m'^3*']'),line=0.75)

##--GSM--##################
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
for(i in 6:1){
  plot(time,latmean(CHL_GSM,loni=loni,lati=LAT[[i]]),type='l',col='dark green',yaxt='n',xlim=timedims,xaxt='n'); 
  axis(side=2,col='dark green'); 
  axis(1,labels=NA,at=ticks)
  
  par(new=TRUE)
  plot(time,latmean(CARBON_GSM,loni=loni,lati=LAT[[i]]),type='l',col='black',xaxt='n',yaxt='n',xlim=timedims)
  axis(4,col='black')
  abline(v=ticks,lty=2)
  mtext(adj=0,labs[i],cex=0.6)
  mtext(paste(round(cor(chl[[i]],carbon[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)
}
axis(1,at=ticks)
mtext(outer=TRUE, side=2, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(outer=TRUE, side=4, expression('Carbon Concentration [mg/m'^3*']'),line=0.75)


###########################################################
## LAG PLOTS ##############################################
###########################################################

labs <- c(expression('00-10'*degree*'N'),expression('10-20'*degree*'N'),expression('20-30'*degree*'N'),
          expression('30-40'*degree*'N'),expression('40-50'*degree*'N'),expression('50-60'*degree*'N'))
yrs <- 2003:2019

pdf('~plots/lag_chl_c.pdf',height=4,width=6)
  par(mfrow=c(1,1),cex.axis=0.8)
  plot(-999,xlim=c(0.5,6.5),ylim=c(-300,100),bty='n',xaxt='n',xlab='',ylab='')
  plotlags(findlag(CHL,CARBON))
  abline(h=0,lty=2)
  axis(side=1,labels=labs,at=1:6)
  mtext(side=2,line=2.5,'Lag of Chlorophyll Maximum [Days]')
dev.off()

pdf('~plots/lag_chl_c_GIOP.pdf',height=4,width=6)
  par(mfrow=c(1,1),cex.axis=0.8)
  plot(-999,xlim=c(0.5,6.5),ylim=c(-300,100),bty='n',xaxt='n',xlab='',ylab='')
  plotlags(findlag(CHL_GIOP,CARBON_GIOP))
  abline(h=0,lty=2)
  axis(side=1,labels=labs,at=1:6)
  mtext(side=2,line=2.5,'Lag of Chlorophyll Maximum [Days]')
dev.off()

pdf('~plots/lag_chl_c_GSM.pdf',height=4,width=6)
  par(mfrow=c(1,1),cex.axis=0.8)
  plot(-999,xlim=c(0.5,6.5),ylim=c(-300,100),bty='n',xaxt='n',xlab='',ylab='')
  plotlags(findlag(CHL_GSM,CARBON_GSM))
  abline(h=0,lty=2)
  axis(side=1,labels=labs,at=1:6)
  mtext(side=2,line=2.5,'Lag of Chlorophyll Maximum [Days]')
dev.off()

###########################################################
## CHL:C vs. MLD-averaged IRRADIANCE ######################
###########################################################
cols <-viridis(12)
labs <- c('00-10 Deg. N', '10-20 Deg. N', '20-30 Deg. N', '30-40 Deg. N', '40-50 Deg. N', '50-60 Deg. N')

pdf('~/dropbox/working/gradients/tzcf/plots/CHLC_vs_E.pdf',height=6,width=5)
par(mfrow=c(3,2),mar=c(2,1,1,4),oma=c(3,3,2,3),xpd=FALSE)
for(i in 6:1){ 
  plot(latmean(I,loni=loni,lati=LAT[[i]]),latmeanratio(CHL,CARBON,loni=loni,lati=LAT[[i]]),pch=19,cex=0.5,col=cols[month(date)],xlim=c(0,300),ylim=c(0,0.05))
  if(i==5){image.plot(matrix(1:12),col=cols,legend.only=TRUE)}
  mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE,expression('Mixed Layer Averaged Irradiance ['*mu*'E/m'^2*'/s]'),line=0.8) 
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 
dev.off()

##--GIOP--#############################
par(mfrow=c(3,2),mar=c(2,1,1,4),oma=c(3,3,2,3),xpd=FALSE)
for(i in 6:1){ 
  plot(latmean(I,loni=loni,lati=LAT[[i]]),latmeanratio(CHL_GIOP,CARBON_GIOP,loni=loni,lati=LAT[[i]]),pch=19,cex=0.5,col=cols[month(date)],xlim=c(0,300),ylim=c(0,0.05))
  if(i==5){image.plot(matrix(1:12),col=cols,legend.only=TRUE)}
  mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE,expression('Mixed Layer Averaged Irradiance ['*mu*'E/m'^2*'/s]'),line=0.8) 
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 

##--GSM--#############################
par(mfrow=c(3,2),mar=c(2,1,1,4),oma=c(3,3,2,3),xpd=FALSE)
for(i in 6:1){ 
  plot(latmean(I,loni=loni,lati=LAT[[i]]),latmeanratio(CHL_GSM,CARBON_GSM,loni=loni,lati=LAT[[i]]),pch=19,cex=0.5,col=cols[month(date)],xlim=c(0,300),ylim=c(0,0.05))
  if(i==5){image.plot(matrix(1:12),col=cols,legend.only=TRUE)}
  mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE,expression('Mixed Layer Averaged Irradiance ['*mu*'E/m'^2*'/s]'),line=0.8) 
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 

##############################################################
## SCATTERPLOTES WITH MLD ####################################
##############################################################
mldin <- seq(0,150,0.1)

labs <- c('a) 50-60 Deg.N','g)',
          'b) 40-50 Deg.N','h)',
          'c) 30-40 Deg.N','i)',
          'd) 20-30 Deg.N','j)',
          'e) 10-20 Deg.N','k)',
          'f) 00-10 Deg.N','l)')

##--GIOP--###########################
pdf('~plots/scatter_chl_c_mld_giop.pdf',height=6,width=6)
par(mfrow=c(6,2),oma=c(3,2,2,2),mar=c(1,4,0,0),cex.axis=0.8,xpd=FALSE)
k <- 1
for(i in 6:1){
  mldi <- latmean(MLD,loni=loni,lati=LAT[[i]])
  ci   <- latmean(CARBON_GIOP,loni=loni,lati=LAT[[i]])
	plot(mldi,ci,pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,50),xaxt='n',xlab='',ylab='')
		linep(mldi,ci,mldin)
		if(i==1) axis(1)
		if(i==3) mtext(side=2,expression('Carbon Concentration [mg/m'^3*']'),line=2,adj=0.335)
		mtext(adj=0,cex=0.6,labs[k])
		k=k+1
		plot(mldi,chli,pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,0.5),xaxt='n',xlab='',ylab='')
		linep(mldi,chli,mldin)
		if(i==1) axis(1)
		if(i==3) mtext(side=2,expression('Chlorophyll Concentration [mg/m'^3*']'),line=2,adj=0.335)
		mtext(adj=0,cex=0.6,labs[k])
		k=k+1
}
mtext(side=1,outer=TRUE,'Mixed Layer Depth [m]',line=1.5)
dev.off()

##--GSM--############################
pdf('~plots/scatter_chl_c_mld_gsm.pdf',height=6,width=6)
par(mfrow=c(6,2),oma=c(3,2,2,2),mar=c(1,4,0,0),cex.axis=0.8,xpd=FALSE)
k <- 1
for(i in 6:1){
  mldi <- latmean(MLD,loni=loni,lati=LAT[[i]])
  ci   <- latmean(CARBON_GSM,loni=loni,lati=LAT[[i]])
  plot(mldi,ci,pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,50),xaxt='n',xlab='',ylab='')
  linep(mldi,ci,mldin)
  if(i==1) axis(1)
  if(i==3) mtext(side=2,expression('Carbon Concentration [mg/m'^3*']'),line=2,adj=0.335)
  mtext(adj=0,cex=0.6,labs[k])
  k=k+1
  plot(mldi,chli,pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,0.5),xaxt='n',xlab='',ylab='')
  linep(mldi,chli,mldin)
  if(i==1) axis(1)
  if(i==3) mtext(side=2,expression('Chlorophyll Concentration [mg/m'^3*']'),line=2,adj=0.335)
  mtext(adj=0,cex=0.6,labs[k])
  k=k+1
}
mtext(side=1,outer=TRUE,'Mixed Layer Depth [m]',line=1.5)
dev.off()


