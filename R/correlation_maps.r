library(reshape)
library(maps)

source("r/load_data.r")

####set box and date parameters####
lon_start  <- -180
lon_end    <- -115
lat_start  <- 0
lat_end    <- 65
date_start <- 2003
date_end   <- 2020

lon <- seq(-180,-115,length.out=130)
lat <- seq(0,65,length.out=130)

ii <- which(time>=(date_start) & time<date_end)

subset_ML     <- MLD[,,ii]
subset_chl    <- CHL[,,ii]
subset_carbon <- CARBON[,,ii]
subset_cbpm   <- CBPM[,,ii]
subset_I      <- I[,,ii]

chl_cor      <- array(numeric(),c(130,130))
carbon_cor   <- array(numeric(),c(130,130))
mu_cor       <- array(numeric(),c(130,130))
cbpm_cor     <- array(numeric(),c(130,130))
I_cor_chl    <- array(numeric(),c(130,130))
I_cor_carbon <- array(numeric(),c(130,130))


#go grid point through grid point and establish correlation, avg derivative value, %missing chl values, and pvalue#

for (p in (1:130)) { print(p)
  for (j in 1:130) {
    if (sum(!is.na(subset_ML[j,p,]))>5 & sum(!is.na(subset_chl[j,p,]))>5 & sum(!is.na(subset_carbon[j,p,]))>5 & sum(!is.na(subset_I[j,p,]))>5) {
      
    chl_cor[j,p]      <- cor(subset_ML[j,p,],subset_chl[j,p,],method = "pearson", use = "pairwise.complete.obs")
    carbon_cor[j,p]   <- cor(subset_ML[j,p,],subset_carbon[j,p,],method = "pearson", use = "pairwise.complete.obs")
    cbpm_cor[j,p]      <- cor(subset_ML[j,p,],subset_cbpm[j,p,],method = "pearson", use = "pairwise.complete.obs")
    I_cor_chl[j,p]    <- cor(subset_I[j,p,],subset_chl[j,p,],method="pearson",use="pairwise.complete.obs")
    I_cor_carbon[j,p] <- cor(subset_I[j,p,],subset_carbon[j,p,],method="pearson",use="pairwise.complete.obs")
    #SODA_interp <- approx(mld_time[MLDi],subset_ML[j,p,],xout = chl_time[chli])$y
      
    # chl_cor[j,p]<-cor(SODA_interp,subset_chl[j,p,],method = "pearson", use = "complete.obs")
    #carbon_cor[j,p]<-cor(SODA_interp,subset_carbon[j,p,],method = "pearson", use = "complete.obs")
    #mu_cor[j,p]<-cor(SODA_interp,subset_mu[j,p,],method = "pearson", use = "complete.obs")
    #NPP_cor[j,p]<-cor(SODA_interp,subset_NPP[j,p,],method = "pearson", use = "complete.obs")
    }
      
    
  }
}


zonal_chl_mld    <- apply(chl_cor,2,mean,na.rm = TRUE)
zonal_carbon_mld <- apply(carbon_cor,2,mean,na.rm = TRUE)
zonal_chl_I      <- apply(I_cor_chl,2,mean,na.rm = TRUE)
zonal_carbon_I   <- apply(I_cor_carbon,2,mean,na.rm = TRUE)


#plot arrays onto map
pdf('plots/correlation_maps.pdf',height=10,width=10)
par(mfrow=c(2,2),mar=c(2,2,2,5),oma=c(2,2,2,2))
image(x=lon,y=lat,carbon_cor,zlim=c(-0.8,0.8),axes=FALSE,col=brewer.pal(n=9,name='RdBu'),xlab='',ylab='')
  map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey',add=TRUE)
  box()
  axis(side=1); axis(side=2); mtext(side=c(3),c('a) MLD vs. Carbon'),adj=0)
  abline(v=seq(-180,-115,10),lty=3); abline(h=seq(0,60,10),lty=3)

image(x=lon,y=lat,chl_cor,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),xlab='',ylab='')
  map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey',add=TRUE)
  box()
  axis(side=1); axis(side=2); mtext(side=c(3),c('b) MLD vs. Chlorophyll'),adj=0)
  image.plot(chl_cor,zlim=c(-.8,.8),legend.only=TRUE,col=brewer.pal(n=11,name='RdBu'))
  abline(v=seq(-180,-115,10),lty=3); abline(h=seq(10,60,10),lty=3)

image(x=lon,y=lat,I_cor_chl,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),xlab='',ylab='')
  map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey',add=TRUE)
  box()
  axis(side=1); axis(side=2); mtext(side=c(3),c('b) MLD-averaged Irradiance vs. Chlorophyll'),adj=0)
  abline(v=seq(-170,-120,10),lty=3); abline(h=seq(10,60,10),lty=3)

image(x=lon,y=lat,I_cor_carbon,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),xlab='',ylab='')
  map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey',add=TRUE)
  box()
  axis(side=1); axis(side=2); mtext(side=c(3),c('a) MLD-averaged Irradiance vs. Carbon'),adj=0)
  abline(v=seq(-170,-120,10),lty=3)
  abline(h=seq(10,60,10),lty=3)
dev.off()



pdf('plots/zonal_correlations.pdf',height=6,width=8)
par(mfrow=c(1,1),oma=c(2,2,2,2))
plot(lat,zonal_chl_mld,'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitute',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonal_carbon_mld, col='darkgreen')
lines(lat,zonal_chl_I, col='green',lty=2)
lines(lat,zonal_carbon_I, col='dark green', lty=2)
  abline(h=0,lty=2)
  legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. Irradiance','Carbon vs. Irradiance'),
         col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2))
dev.off()
