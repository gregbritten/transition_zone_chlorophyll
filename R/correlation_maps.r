rm(list=ls())
library(reshape)
library(maps)

source("r/load_data.r")
source("r/functions.r")

#plot arrays onto map
pdf('plots/correlation_maps.pdf',height=10,width=10)
par(mfrow=c(2,2),mar=c(2,2,2,5),oma=c(2,2,2,2))
plotcormap(cormap(CARBON,MLD),lab='a) MLD vs. Carbon')
plotcormap(cormap(CHL,MLD),lab='b) MLD vs. Chlorophyll')
plotcormap(cormap(CARBON,I),lab='c) MLD-averaged Irradiance vs. Carbon')
plotcormap(cormap(CHL,I),lab='c) MLD-averaged Irradiance vs. Chlorophyll')
dev.off()

par(mfrow=c(2,2),mar=c(2,2,2,5),oma=c(2,2,2,2))
plotcormap(cormap(CARBON_GIOP,MLD),lab='a) MLD vs. Carbon')
plotcormap(cormap(CHL_GIOP,MLD),lab='b) MLD vs. Chlorophyll')
plotcormap(cormap(CARBON_GIOP,I),lab='c) MLD-averaged Irradiance vs. Carbon')
plotcormap(cormap(CHL_GIOP,I),lab='c) MLD-averaged Irradiance vs. Chlorophyll')

par(mfrow=c(2,2),mar=c(2,2,2,5),oma=c(2,2,2,2))
plotcormap(cormap(CARBON_GSM,MLD),lab='a) MLD vs. Carbon')
plotcormap(cormap(CHL_GSM,MLD),lab='b) MLD vs. Chlorophyll')
plotcormap(cormap(CARBON_GSM,I),lab='c) MLD-averaged Irradiance vs. Carbon')
plotcormap(cormap(CHL_GSM,I),lab='c) MLD-averaged Irradiance vs. Chlorophyll')




pdf('plots/zonal_correlations.pdf',height=6,width=8)
par(mfrow=c(1,1),oma=c(2,2,2,2))
plot(lat,zonalcor(cormap(CHL,MLD)),'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitute',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonalcor(cormap(CARBON,MLD)), col='darkgreen')
lines(lat,zonalcor(cormap(CHL,I)), col='green',lty=2)
lines(lat,zonalcor(cormap(CARBON,I)), col='dark green', lty=2)
  abline(h=0,lty=2)
  legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. Irradiance','Carbon vs. Irradiance'),
         col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2))
dev.off()


par(mfrow=c(1,1),oma=c(2,2,2,2))
plot(lat,zonalcor(cormap(CHL_GIOP,MLD)),'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitute',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonalcor(cormap(CARBON_GIOP,MLD)), col='darkgreen')
lines(lat,zonalcor(cormap(CHL_GIOP,I)), col='green',lty=2)
lines(lat,zonalcor(cormap(CARBON_GIOP,I)), col='dark green', lty=2)
abline(h=0,lty=2)
legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. Irradiance','Carbon vs. Irradiance'),
       col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2))


par(mfrow=c(1,1),oma=c(2,2,2,2))
plot(lat,zonalcor(cormap(CHL_GSM,MLD)),'l',ylim = c(-.8,.8),xlim=c(0,70),col='green',xlab = 'Latitute',ylab = 'Zonal Average Correlations',bty='n',yaxt='n')
axis(side=2,at=seq(-0.8,0.8,0.2))
lines(lat,zonalcor(cormap(CARBON_GSM,MLD)), col='darkgreen')
lines(lat,zonalcor(cormap(CHL_GSM,I)), col='green',lty=2)
lines(lat,zonalcor(cormap(CARBON_GSM,I)), col='dark green', lty=2)
abline(h=0,lty=2)
legend("topright",c('Chlorophyll vs. MLD','Carbon vs. MLD','Chlorophyll vs. Irradiance','Carbon vs. Irradiance'),
       col=c('green','darkgreen','green','dark green'), bty='n', lty=c(1,1,2,2))




