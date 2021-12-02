
rm(list=ls())
library(viridis)
library(fields)

source("~/dropbox/working/gradients/tzcf/github_mac/greg/load_data.r")
mnth_ocn <- month(ocnprd_date)
#mnth_mld <- month(mld_date)

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

nlat  <- length(LAT)
nmnth <- 12

chl=mld=c=chlc=mu=npp=par=kd=e=dn=parmld <- matrix(NA,nrow=nlat,ncol=nmnth)
for(i in 1:nlat){
  for(k in 1:nmnth){
    
    tmpchl                <- CHL[loni,LAT[[i]],mnth_ocn==k]
    chl[i,k]              <- mean(tmpchl,na.rm=TRUE)
    
    tmpmld                <- MLD[loni,LAT[[i]],mnth_ocn==k]
    mld[i,k]              <- mean(tmpmld,na.rm=TRUE)
    
    tmpc                  <- C[loni,LAT[[i]],mnth_ocn==k]
    c[i,k]                <- mean(tmpc,na.rm=TRUE)
    
    chlc[i,k]             <- chl[i,k]/c[i,k] 
    
    tmpmu                 <- MU[loni,LAT[[i]],mnth_ocn==k]
    mu[i,k]               <- mean(tmpmu, na.rm=TRUE)
    
    npp[i,k]              <- mu[i,k]*c[i,k]
    
    tmppar                <- PAR[loni,LAT[[i]],mnth_ocn==k]
    par[i,k]              <- mean(tmppar, na.rm=TRUE)
    
    tmpk                  <- K[loni,LAT[[i]],mnth_ocn==k]
    kd[i,k]               <- mean(tmpk, na.rm=TRUE)
    
    tmpparmld             <- I[loni,LAT[[i]],mnth_ocn==k]
    parmld[i,k]           <- mean(tmpparmld,na.rm=TRUE)
    
    e[i,k]    <- (1/(kd[i,k]*mld[i,k]))*par[i,k]*exp(-(kd[i,k]*mld[i,k]))*(1-exp(-(kd[i,k]*mld[i,k])))
  }
}


pdf('~/dropbox/working/gradients/tzcf/plots/surfI_MLD_MLDI_lats.pdf',height=6,width=7)
par(mfrow=c(3,1),mar=c(2,4,1,8),oma=c(2,2,2,6)) 
  matplot(t(par),type='l',col=viridis(6),lty=1,ylab='',xlab='')
    mtext('Surface Irradiance',side=2,line=3.5)
    mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
    image.plot(matrix(c(10,60)), legend.only=TRUE,col=viridis(6))  
  matplot(t(mld),type='l',col=viridis(6),lty=1,ylab='',xlab='')
    mtext('MLD',side=2,line=3.5)
    mtext('[m]',side=2,line=2.5,cex=0.7)
  matplot(t(e),type='l',col=viridis(6),lty=1,ylab='',xlab='')
    mtext('MLD-average Irradiance',side=2,line=3.5)
    mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
dev.off()




