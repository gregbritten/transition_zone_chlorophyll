

rm(list=ls())
library(viridis)
library(fields)

#source('d:/dropbox/working/gradients/tzcf/github/greg/load_data.r')
source("~/dropbox/working/gradients/tzcf/github_mac/greg/load_data.r")

lon <- seq(-180,-115,length.out=130)
lat <- seq(0,65,length.out=130)

loni <- which(lon > -180 & lon < -115)
LAT <- list()
LAT[['lat00_10']] <- which(lat >=0  & lat <10)
LAT[['lat10_20']] <- which(lat >=10 & lat <20)
LAT[['lat20_30']] <- which(lat >=20 & lat <30)
LAT[['lat30_40']] <- which(lat >=30 & lat <40)
LAT[['lat40_50']] <- which(lat >=40 & lat <50)
LAT[['lat50_60']] <- which(lat >=50 & lat <60)

chl=mld=c=chlc=mu=npp=par=k=e <- list()
for(i in 1:6){
  tmpchl                <- CHL[loni,LAT[[i]],]
  chl[[names(LAT)[i]]]  <- apply(tmpchl,3,function(x) mean(x,na.rm=TRUE))
  
  tmpmld                <- MLD[loni,LAT[[i]],]
  mld[[names(LAT)[i]]]  <- apply(tmpmld,3,function(x) mean(x,na.rm=TRUE))
  
  tmpc                  <- C[loni,LAT[[i]],]
  c[[names(LAT)[i]]]    <- apply(tmpc,3,function(x) mean(x, na.rm=TRUE))
  
  chlc[[names(LAT)[i]]] <- chl[[names(LAT)[i]]]/c[[names(LAT)[i]]] 
  
  tmpmu                 <- MU[loni,LAT[[i]],]
  mu[[names(LAT)[i]]]   <- apply(tmpmu,3,function(x) mean(x, na.rm=TRUE))
  
  npp[[names(LAT)[i]]]  <- mu[[names(LAT)[i]]]*c[[names(LAT)[i]]]
  
  tmppar                <- PAR[loni,LAT[[i]],]
  par[[names(LAT)[i]]]  <- apply(tmppar,3,function(x) mean(x, na.rm=TRUE))
  
  tmpk                  <- K[loni,LAT[[i]],]
  k[[names(LAT)[i]]]    <- apply(tmpk,3,function(x) mean(x, na.rm=TRUE))
  
  e[[names(LAT)[i]]]    <- (1/(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]]))*par[[names(LAT)[i]]]*exp(-(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]]))*(1-exp(-(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]])))
}

###########################################################
## CHL:C vs. MLD-averaged IRRADIANCE ######################
###########################################################
cols <-viridis(12)
labs <- c('00-10 Deg. N', '10-20 Deg. N', '20-30 Deg. N', '30-40 Deg. N', '40-50 Deg. N', '50-60 Deg. N')

pdf('~/dropbox/working/gradients/tzcf/plots/CHLC_vs_E.pdf',height=6,width=5)
par(mfrow=c(3,2),mar=c(2,1,1,4),oma=c(3,3,2,3),xpd=FALSE)
for(i in 6:1){ 
  plot(e[[i]],chlc[[i]],pch=19,cex=0.5,col=cols[month(ocnprd_date)],xlim=c(0,300),ylim=c(0,0.05))
  if(i==5){image.plot(matrix(1:12),col=cols,legend.only=TRUE)}
  mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE,expression('Mixed Layer Averaged Irradiance ['*mu*'E/m'^2*'/s]'),line=0.8) 
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 
dev.off()
