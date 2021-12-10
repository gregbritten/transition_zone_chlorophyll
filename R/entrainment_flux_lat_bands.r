
rm(list=ls())
library(viridis)
library(fields)

source("r/load_data.r")
source("r/functions.r")

mnth <- month(date)

e    <- (1/(monthmean(KD)*monthmean(MLD)))*monthmean(PAR)*exp(-(monthmean(KD)*monthmean(MLD)))*(1-exp(-(monthmean(KD)*monthmean(MLD))))

##--COMPUTE N FLUXES--###############################
zz  <- 1:30
zs  <- woa_z[zz]
mld <- monthmean(MLD) 

dn <- matrix(NA,nrow=6,ncol=12)
for(i in 1:6){
  for(k in 1:12){
      tmpN    <- apply(N[,LAT[[i]],zz,k],3,function(x) mean(x,na.rm=TRUE))
      f       <- approxfun(x=zs,y=tmpN)
      dn[i,k] <- f(mld[i,k]) - mean(f(1:mld[i,k]),na.rm=TRUE)
  }
}
  
dmld <- apply(cbind(mld[,12],mld),1,diff)
nflux <- t(dn)*dmld  
nflux[nflux<0] <- 0

######################################################
## PLOT ##############################################
######################################################
labs <- c('00-10 Deg. N', '10-20 Deg. N', '20-30 Deg. N', '30-40 Deg. N', '40-50 Deg. N', '50-60 Deg. N')

pdf('~/dropbox/working/gradients/tzcf/plots/mld_dmlddt_N0_FN_lats_03_29_2021.pdf',height=7,width=7)
par(mfrow=c(4,1),mar=c(2,4,1,8),oma=c(2,2,2,6)) 
matplot(t(monthmean(MLD)),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
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
matplot(t(monthmean(PAR)),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('Surface Irradiance',side=2,line=3.5)
  mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
image.plot(matrix(c(10,60)), legend.only=TRUE,col=viridis(6))  
matplot(t(monthmean(MLD)),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('MLD',side=2,line=3.5)
  mtext('[m]',side=2,line=2.5,cex=0.7)
matplot(t(e),type='l',col=viridis(6),lty=c(1,1,1,2,1,1),ylab='',xlab='')
  mtext('MLD-average Irradiance',side=2,line=3.5)
  mtext(expression(mu*'E/m'^2*'/s'),side=2,line=2.25,cex=0.7)
dev.off()

