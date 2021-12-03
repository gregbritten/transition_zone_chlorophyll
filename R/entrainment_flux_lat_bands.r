
rm(list=ls())
library(viridis)
library(fields)

source("r/load_data.r")

load(file="~/dropbox/working/gradients/tzcf/data/WOA_N.rdata")
load(file="~/dropbox/working/gradients/tzcf/data/WOA_P.rdata")
load(file="~/dropbox/working/gradients/tzcf/data/WOA_depth.rdata")

lat  <- seq(0,65,length.out=130)
lon  <- seq(-180,-115,length.out=130)
mnth <- month(date)

LAT <- list()
LAT[['lat00_10']] <- which(lat >=0  & lat <10)
LAT[['lat10_20']] <- which(lat >=10 & lat <20)
LAT[['lat20_30']] <- which(lat >=20 & lat <30)
LAT[['lat30_40']] <- which(lat >=30 & lat <40)
LAT[['lat40_50']] <- which(lat >=40 & lat <50)
LAT[['lat50_60']] <- which(lat >=50 & lat <60)

zz <- 1:30
zs <- z[zz]

chl=mld=zmld=c=chlc=mu=npp=par=kd=e=dn=parmld <- matrix(NA,nrow=6,ncol=12)
for(i in 1:6){
  for(k in 1:12){
  
      tmpchl                <- CHL[,LAT[[i]],mnth==k]
      chl[i,k]              <- mean(tmpchl,na.rm=TRUE)
      
      tmpmld                <- MLD[,LAT[[i]],mnth==k]
      mld[i,k]              <- mean(tmpmld,na.rm=TRUE)
      
      tmpc                  <- CARBON[,LAT[[i]],mnth==k]
      c[i,k]                <- mean(tmpc,na.rm=TRUE)
      
      chlc[i,k]             <- chl[i,k]/c[i,k] 
      
      #npp[i,k]              <- mu[i,k]*c[i,k]
      
      tmppar                <- PAR[,LAT[[i]],mnth==k]
      par[i,k]              <- mean(tmppar, na.rm=TRUE)
      
      tmpk                  <- KD[,LAT[[i]],mnth==k]
      kd[i,k]               <- mean(tmpk, na.rm=TRUE)
      
      tmpparmld             <- I[,LAT[[i]],mnth==k]
      parmld[i,k]           <- mean(tmpparmld,na.rm=TRUE)
      
      e[i,k]    <- (1/(kd[i,k]*mld[i,k]))*par[i,k]*exp(-(kd[i,k]*mld[i,k]))*(1-exp(-(kd[i,k]*mld[i,k])))
      
      
      tmpN    <- apply(N[,LAT[[i]],zz,k],3,function(x) mean(x,na.rm=TRUE))
      f       <- approxfun(x=zs,y=tmpN)
      
      zmld[i,k] <- mld[i,k]
      dn[i,k]   <- f(mld[i,k]) - mean(f(1:mld[i,k]),na.rm=TRUE)
  }
}
  
dmld <- apply(cbind(zmld[,12],zmld),1,diff)
nflux <- t(dn)*dmld  
nflux[nflux<0] <- 0
nflux2 <- t(dn[,c(12,1:11)])*dmld  
nflux2[nflux2<0] <- 0
nflux3 <- t(dn[,c(2:12,1)])*dmld  
nflux3[nflux3<0] <- 0

par(mfrow=c(2,1),mar=c(2,4,1,8),oma=c(2,2,2,6)) 
matplot(nflux,type='l',ylim=c(0,60),col=viridis(6),lty=1,ylab='',xlab='')  
  mtext('Nitrate gradient from same month')
  image.plot(matrix(c(10,60)), legend.only=TRUE,col=viridis(6))  
matplot(nflux2,type='l',ylim=c(0,60),col=viridis(6),lty=1,ylab='',xlab='')  
  mtext('Nitrate gradient from previous month')
#matplot(nflux3,type='l',ylim=c(0,60),col=viridis(6),lty=1,ylab='',xlab='')  
  mtext(expression(italic('F'['N'])),side=2,line=-1,outer=TRUE)
  mtext(expression('[mmol/m'^2*'/month]'),side=2,line=-2,cex=0.7,outer=TRUE)
 
labs <- c('00-10 Deg. N', '10-20 Deg. N', '20-30 Deg. N', '30-40 Deg. N', '40-50 Deg. N', '50-60 Deg. N')
  

pdf('~/dropbox/working/gradients/tzcf/plots/mld_dmlddt_N0_FN_lats_03_29_2021.pdf',height=7,width=7)
par(mfrow=c(4,1),mar=c(2,4,1,8),oma=c(2,2,2,6)) 
matplot(t(zmld),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('MLD',side=2,line=3.5)
  mtext('[m]',side=2,line=2.5,cex=0.7)
  image.plot(matrix(c(10,60)), legend.only=TRUE,col=viridis(6))  
matplot(dmld,type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('dMLD/dt',side=2,line=3.5)
  mtext('[m/month]',side=2,line=2.5,cex=0.7)
  
abline(h=0,lty=1,lwd=1)
matplot(t(dn),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext(expression(italic('N'['0']~'- N')),side=2,line=3.5)
  mtext(expression('[mmol/m'^3*']'),side=2,line=2.25,cex=0.7)
matplot(nflux,type='l',ylim=c(0,60),col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')  
  mtext(expression(italic('F'['N'])),side=2,line=3.5)
  mtext(expression('[mmol/m'^2*'/month]'),side=2,line=2.25,cex=0.7)
dev.off()
  

pdf('~/dropbox/working/gradients/tzcf/plots/surfI_MLD_MLDI_lats_03_29_2021.pdf',height=6,width=7)
par(mfrow=c(3,1),mar=c(2,4,1,8),oma=c(2,2,2,6)) 
matplot(t(par),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('Surface Irradiance',side=2,line=3.5)
  mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
image.plot(matrix(c(10,60)), legend.only=TRUE,col=viridis(6))  
matplot(t(zmld),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('MLD',side=2,line=3.5)
  mtext('[m]',side=2,line=2.5,cex=0.7)
matplot(t(e),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('MLD-average Irradiance',side=2,line=3.5)
  mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
#matplot(t(parmld),type='l',col=viridis(6),lty=1)
dev.off()




 
write.csv(file='~/dropbox/working/gradients/tzcf/data/N_FLUX.csv',nflux)  





for(k in 1:12){
  tmpN    <- colMeans(N[i,j,zz,k],na.rm=TRUE)
  mlds <- MLD[i,j,mnth_mld==k]
  mld  <- mean(mlds,na.rm=TRUE)
  
  f <- approxfun(x=zs,y=tmpN)
  
  ZMLD[k,l] <- mld
  DN[k,l]   <- f(mld) - mean(f(1:mld),na.rm=TRUE)
  
  #D[[l]][[k]] <- list()
  #D[[l]][[k]]$N    <- colMeans(N[i,j,zz,k])
  #D[[l]][[k]]$mlds     <- MLD[i,j,mnth_mld==k]
  #round(lats[j],1)
  pars <- PAR[i,j,mnth[1:812]==k]
  pars <- pars/86400*1E6/4.6
  PARS[k,l] <- mean(pars,na.rm=TRUE)
  ks   <- K490[i,j,mnth[1:812]==k]
  K[k,l] <- mean(ks,na.rm=TRUE)
  
  E[k,l]    <- (1/(K[k,l]*mld))*PARS[k,l]*exp(-(K[k,l]*mld))*(1-exp(-(K[k,l]*mld)))
  
  #ks_mean <- 0.1
  #light <- par_mean*exp(-ks_mean*z[zz])
  #zcriti <- which((light - 21)^2 == min((light - 21)^2))
}


