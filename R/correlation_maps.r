### This script calculates the correlation between SODA data and chl readings. 
####Load Libraries ####

library(lubridate)
library(dplyr)
library(tidyr)
library(reshape)
library(sp) 
library(ggplot2)
library(fields)
library(viridis)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(zoo)
setwd('~/dropbox/working/gradients/tzcf/data/')

####import SODA data####
load('SODA_lats.rdata')
load('SODA_lons.rdata')
load('SODA_area.rdata')
#load('SODA_time.rdata')
#mld_time<-time
load('OCNPRD_time.rdata')
load('OCNPRD_MLD_time.rdata')
mld_time <- decimal_date(ymd(mld_time))

chl_time <- decimal_date(ymd(ocnprd_date))
#mld_time <- chl_time

#load('ML.rdata')
load('OCNPRD_MLD.rdata')
ML<- DAT
#load('PAR_MLD.rdata')
#ML<-I

load('OCNPRD_growth_rate.rdata')
mu <- DAT

load('OCNPRD_CHL_gapfilled.rdata')
chl <- DAT

load('OCNPRD_carbon.rdata')
carbon <- DAT

load('PAR_MLD.rdata')

NPP <- mu*carbon




####set box and date parameters####
lon_start<- -180
lon_end <- -115
lat_start <- 0
lat_end <- 60
date_start <- 2003
date_end <- 2019

MLDi <- which(mld_time>=(date_start)&mld_time<date_end)
chli<-which(chl_time>=(date_start)&chl_time<(date_end))

subset_ML <- ML[,,MLDi]
subset_chl <- chl[,,chli]
subset_carbon <- carbon[,,chli]
subset_mu <- mu[,,chli]
subset_NPP <- NPP[,,chli]
subset_I <- I[,,chli]


chl_cor <- array(numeric(),c(130,130))
carbon_cor <- array(numeric(),c(130,130))
mu_cor <- array(numeric(),c(130,130))
NPP_cor <- array(numeric(),c(130,130))
I_cor_chl <- array(numeric(),c(130,130))
I_cor_carbon <- array(numeric(),c(130,130))


#go grid point through grid point and establish correlation, avg derivative value, %missing chl values, and pvalue#

for (p in (1:130)) { print(p)
  for (j in 1:130) {
    

    if (length(which(!is.na(subset_ML[j,p,])))>5&length(which(!is.na(subset_chl[j,p,])))>5&length(which(!is.na(subset_carbon[j,p,])))>5) {
      
      chl_cor[j,p]<-cor(subset_ML[j,p,],subset_chl[j,p,],method = "pearson", use = "complete.obs")
    carbon_cor[j,p]<-cor(subset_ML[j,p,],subset_carbon[j,p,],method = "pearson", use = "complete.obs")
    mu_cor[j,p]<-cor(subset_ML[j,p,],subset_mu[j,p,],method = "pearson", use = "complete.obs")
    NPP_cor[j,p]<-cor(subset_ML[j,p,],subset_NPP[j,p,],method = "pearson", use = "complete.obs")
    I_cor_chl[j,p] <- cor(subset_I[j,p,],subset_chl[j,p,],method="pearson",use='complete.obs')
    I_cor_carbon[j,p] <- cor(subset_I[j,p,],subset_carbon[j,p,],method="pearson",use='complete.obs')
    #SODA_interp <- approx(mld_time[MLDi],subset_ML[j,p,],xout = chl_time[chli])$y
      
    
      # chl_cor[j,p]<-cor(SODA_interp,subset_chl[j,p,],method = "pearson", use = "complete.obs")
    #carbon_cor[j,p]<-cor(SODA_interp,subset_carbon[j,p,],method = "pearson", use = "complete.obs")
    #mu_cor[j,p]<-cor(SODA_interp,subset_mu[j,p,],method = "pearson", use = "complete.obs")
    #NPP_cor[j,p]<-cor(SODA_interp,subset_NPP[j,p,],method = "pearson", use = "complete.obs")
    }
      
    
  }
}


chl_avg<- apply(chl_cor,2,mean,na.rm = TRUE)
carbon_avg <- apply(carbon_cor,2,mean,na.rm = TRUE)
mu_avg <- apply(mu_cor,2,mean,na.rm = TRUE)
NPP_avg<-apply(NPP_cor,2,mean,na.rm = TRUE)


#plot arrays onto map
pdf('~/dropbox/working/gradients/tzcf/plots/carbon_cor.pdf',8)
#par(mar=c(2,2,2,7))
#layout(matrix(c(1,2),ncol = 2),widths=c(3,1))
map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey')
box()
image(x=lon,y=lat,carbon_cor,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),add=TRUE)
axis(side=1); axis(side=2); mtext(side=c(3),c('a) MLD vs. Carbon'),adj=0)
image.plot(carbon_cor,zlim=c(-.8,.8),legend.only=TRUE,col=brewer.pal(n=11,name='RdBu'))
abline(v=seq(-170,-120,10),lty=3)
abline(h=seq(10,60,10),lty=3)
dev.off()

#plot arrays onto map
pdf('~/dropbox/working/gradients/tzcf/plots/chl_cor.pdf',8)

#layout(matrix(c(1,2),ncol = 2),widths=c(3,1))
map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey')
box()
image(x=lon,y=lat,chl_cor,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),add=TRUE)
axis(side=1); axis(side=2); mtext(side=c(3),c('b) MLD vs. Chlorophyll'),adj=0)
image.plot(chl_cor,zlim=c(-.8,.8),legend.only=TRUE,col=brewer.pal(n=11,name='RdBu'))
abline(v=seq(-170,-120,10),lty=3)
abline(h=seq(10,60,10),lty=3)
dev.off()


pdf('~/dropbox/working/gradients/tzcf/plots/I_chl_cor.pdf',8)
map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey')
box()
image(x=lon,y=lat,I_cor_chl,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),add=TRUE)
axis(side=1); axis(side=2); mtext(side=c(3),c('b) MLD-averaged Irradiance vs. Chlorophyll'),adj=0)
image.plot(chl_cor,zlim=c(-.8,.8),legend.only=TRUE,col=brewer.pal(n=11,name='RdBu'))
abline(v=seq(-170,-120,10),lty=3)
abline(h=seq(10,60,10),lty=3)
dev.off()


pdf('~/dropbox/working/gradients/tzcf/plots/I_cor_carbon.pdf',8)
map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey')
box()
image(x=lon,y=lat,I_cor_carbon,zlim=c(-.8,.8),axes = FALSE,col=brewer.pal(n=9,name='RdBu'),add=TRUE)
axis(side=1); axis(side=2); mtext(side=c(3),c('a) MLD-averaged Irradiance vs. Carbon'),adj=0)
image.plot(chl_cor,zlim=c(-.8,.8),legend.only=TRUE,col=brewer.pal(n=11,name='RdBu'))
abline(v=seq(-170,-120,10),lty=3)
abline(h=seq(10,60,10),lty=3)
dev.off()



pdf('avg_correlations.pdf',8,6)
plot(lat,chl_avg,'l',ylim = c(-.8,.8),col='green',xlab = 'Latitute',ylab = 'Zonal Average Correlations')
lines(lat,carbon_avg,col='darkgreen')
lines(lat,NPP_avg,col='blue')
lines(lat,mu_avg,col = 'lightblue')
abline(h=0,lty=2)
legend("topright",c('Chlorophyll','Carbon','NPP','Growth Rate'),col = c('green','darkgreen','blue','lightblue'),bty='n',lty=1)
  dev.off()
