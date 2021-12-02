
rm(list=ls())
library(viridis)
library(fields)
library(RColorBrewer)
#source('d:/dropbox/working/gradients/tzcf/github/greg/load_data.r')
source("~/dropbox/working/gradients/tzcf/github/greg/load_data.r")

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

chl=mld=c=chlc=mu=npp=par=k=e=vgpm=cafe=vgpmep <- list()
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
	
	#npp[[names(LAT)[i]]]  <- mu[[names(LAT)[i]]]*c[[names(LAT)[i]]]
	
	
  tmpcbpm               <- CBPM[loni,LAT[[i]],]
	npp[[names(LAT)[i]]]  <- apply(tmpcbpm,3,function(x) mean(x,na.rm=TRUE))

	tmpvgpm               <- VGPM[loni,LAT[[i]],]
	vgpm[[names(LAT)[i]]] <- apply(tmpvgpm,3,function(x) mean(x,na.rm=TRUE))
	
	tmpcafe               <- CAFE[loni,LAT[[i]],]
	cafe[[names(LAT)[i]]] <- apply(tmpcafe,3,function(x) mean(x,na.rm=TRUE))

	tmpvgpmep             <- EPPLEY_VGPM[loni,LAT[[i]],]
  vgpmep[[names(LAT[i])]]<- apply(tmpvgpmep,3,function(x) mean(x,na.rm=TRUE))
			
	tmppar                <- PAR[loni,LAT[[i]],]
	par[[names(LAT)[i]]]  <- apply(tmppar,3,function(x) mean(x, na.rm=TRUE))
	
	tmpk                  <- K[loni,LAT[[i]],]
	k[[names(LAT)[i]]]    <- apply(tmpk,3,function(x) mean(x, na.rm=TRUE))
	
	e[[names(LAT)[i]]]    <- (1/(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]]))*par[[names(LAT)[i]]]*exp(-(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]]))*(1-exp(-(k[[names(LAT)[i]]]*mld[[names(LAT)[i]]])))
}

labs <- c('f) 00-10 Deg. N', 'e) 10-20 Deg. N', 'd) 20-30 Deg. N', 'c) 30-40 Deg. N', 'b) 40-50 Deg. N', 'a) 50-60 Deg. N')

crit_depth <- 21*4.57
 
timedims <- c(2003,2021)

ticks <- 2003:2020

cols <- viridis(13)

linep <- function(x,y,xin){
	fit <- lm(y ~ x)
	lines(xin,predict(fit,newdata=list(x=xin)),col='dark grey')
}


eee <- as.numeric(unlist(e))
eee[log10(eee)< -1.5] <- 10^(-1.5)
eee[log10(eee)> 1.5] <- 10^(1.5)
cols <- as.character(cut(log10(eee),length(eee),labels=viridis(length(eee))))

COLS <- list()
n0 <- 1
n1 <- 1
for(i in 1:6){
	n <- length(mu[[i]])
	n0 <- n1+1
	n1 <- n1 + n -1
	COLS[[i]] <- cols[n0:n1]
}
	
##########################################################
## GROWTH RATE VS CHL:C ##################################
##########################################################
muin <- seq(0,1.4,0.01)
pdf('~/dropbox/working/gradients/tzcf/plots/chl_c_mu.pdf',height=6,width=6)
par(mfrow=c(3,2),mar=c(2,2,2,4),oma=c(3,3,2,2),xpd=TRUE,cex.axis=0.8,cex.lab=.8)
n0 <- 1
n1 <- length(eee)
for(i in 6:1){ 
	#cols <- as.character(cut(log10(e[[i]]),length(e[[i]]),labels=viridis(length(e[[i]]))))
	#plot(mu[[i]],chlc[[i]],pch=19,cex=0.5,col=cols[month(ocnprd_date)])
	plot(mu[[i]],chlc[[i]],pch=19,cex=0.5,col=COLS[[i]],xlim=c(0,1.4),ylim=c(0,0.08))
	#plot(mu[[i]],chlc[[i]],pch=19,cex=0.5,col=cols)
	#linep(mu[[i]],chlc[[i]],muin)
	#if(i==2){image.plot(matrix(as.numeric(unlist(e))),col=cols,legend.only=TRUE)}
	image.plot(matrix(c(-1.5,1.5)),col=viridis(12),legend.only=TRUE)
	mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE, expression(mu~'[day'^{-1}*']'),line=1)
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 
dev.off()

#####################################################################
## COMPARISON OF NPP PRODUCTS #######################################
#####################################################################
cols <- brewer.pal(4,"Dark2")
pdf('~/dropbox/working/gradients/tzcf/plots/chl_npp_comparison.pdf',height=10,width=8)
par(mfrow=c(12,1),mar=c(1,2,1,2),oma=c(2,3,2,5),cex.axis=0.7)
plot(chl_time,npp[[6]],type='l',col=cols[1],xaxt='n',xlim=timedims,yaxt='n',ylim=range(c(npp[6],cafe[6],vgpm[6],vgpmep[6]),na.rm=TRUE))
lines(chl_time,vgpm[[6]],type='l',col=cols[2],xaxt='n',xlim=timedims,yaxt='n')
lines(chl_time,cafe[[6]],type='l',col=cols[3],xaxt='n',xlim=timedims,yaxt='n')
lines(chl_time,vgpmep[[6]],type='l',col=cols[4],xaxt='n',xlim=timedims,yaxt='n')
legend('topright',legend=c('CBPM','VGPM','CAFE','VGPM_EP'),col=cols,lty=1,cex=0.6,bty='n')
mtext(paste(round(cor(chl[[6]],npp[[i]],use='pairwise.complete.obs'),3)),adj=0.7,cex=0.8,col=cols[1])	
mtext(paste(round(cor(chl[[6]],vgpm[[i]],use='pairwise.complete.obs'),3)),adj=0.8,cex=0.8,col=cols[2])	
mtext(paste(round(cor(chl[[6]],cafe[[i]],use='pairwise.complete.obs'),3)),adj=0.9,cex=0.8,col=cols[3])	
mtext(paste(round(cor(chl[[6]],vgpmep[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.8,col=cols[4])	

mtext(adj=0,labs[6],cex=0.6)
axis(1,labels=NA,at=ticks)
abline(v=ticks,lty=2)
axis(side=2)
plot(chl_time,chl[[6]],type='l',xlim=timedims,xaxt='n',yaxt='n')
axis(1,labels=NA,at=ticks)
abline(v=ticks,lty=2)
axis(side=4)
for(i in 5:1){
  plot(chl_time,npp[[i]],type='l',col=cols[1],xlim=timedims,xaxt='n',yaxt='n',ylim=range(c(npp[i],vgpm[i],cafe[i],vgpmep[i]))); 
  lines(chl_time,vgpm[[i]],type='l',col=cols[2],xaxt='n',xlim=timedims,yaxt='n')
  lines(chl_time,cafe[[i]],type='l',col=cols[3],xaxt='n',xlim=timedims,yaxt='n')
  lines(chl_time,vgpmep[[i]],type='l',col=cols[4],xaxt='n',xlim=timedims,yaxt='n')
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=2)
  par(new=TRUE)	
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chl[[i]],npp[[i]],use='pairwise.complete.obs'),3)),adj=0.7,cex=0.8,col=cols[1])	
  mtext(paste(round(cor(chl[[i]],vgpm[[i]],use='pairwise.complete.obs'),3)),adj=0.8,cex=0.8,col=cols[2])	
  mtext(paste(round(cor(chl[[i]],cafe[[i]],use='pairwise.complete.obs'),3)),adj=0.9,cex=0.8,col=cols[3])	
  mtext(paste(round(cor(chl[[i]],vgpmep[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.8,col=cols[4])	
plot(chl_time,chl[[i]],type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
axis(1,labels=NA,at=ticks)
abline(v=ticks,lty=2)
axis(side=4)

}
axis(1,at=ticks)
mtext(side=4,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'),line=0.75)
mtext(side=2,outer=TRUE, expression('Net Primary Productivity [mg/m'^2*'/day]'))
dev.off()



#####################################################################
## TIME SERIES PLOTS - CHL VS NPP ###################################
#####################################################################
#pdf('d:/dropbox/working/gradients/tzcf/plots/chl_npp.pdf',height=6,width=5)
pdf('~/dropbox/working/gradients/tzcf/plots/chl_npp.pdf',height=6,width=5)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
plot(chl_time,npp[[6]],type='l',col='black',xaxt='n',xlim=timedims,yaxt='n')
	mtext(adj=0,labs[6],cex=0.6)
	axis(1,labels=NA,at=ticks)
	abline(v=ticks,lty=2)
	axis(side=4)
par(new=TRUE)	
plot(chl_time,chl[[6]],type='l',xaxt='n',xlim=timedims,col='dark green',yaxt='n')
	mtext(paste(round(cor(chl[[6]],npp[[6]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
	axis(side=2,col='dark green')
for(i in 5:1){
	plot(chl_time,npp[[i]],type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
	abline(v=ticks,lty=2)
	axis(1,labels=NA,at=ticks)
	axis(side=4)
par(new=TRUE)	
plot(chl_time,chl[[i]],type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
	mtext(adj=0,labs[i],cex=0.6)
	axis(1,labels=NA,at=ticks)
	mtext(paste(round(cor(chl[[i]],npp[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
	axis(side=2,col='dark green')
}
axis(1,at=ticks)
mtext(side=2,outer=TRUE, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(side=4,outer=TRUE, expression('Net Primary Productivity [mg/m'^2*'/day]'),line=0.75)
dev.off()

pdf('~/dropbox/working/gradients/tzcf/plots/chl_npp_vgpm.pdf',height=6,width=5)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)
plot(chl_time,vgpm[[6]],type='l',col='black',xaxt='n',xlim=timedims,yaxt='n')
mtext(adj=0,labs[6],cex=0.6)
axis(1,labels=NA,at=ticks)
abline(v=ticks,lty=2)
axis(side=4)
par(new=TRUE)	
plot(chl_time,chl[[6]],type='l',xaxt='n',xlim=timedims,col='dark green',yaxt='n')
mtext(paste(round(cor(chl[[6]],npp[[6]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
axis(side=2,col='dark green')
for(i in 5:1){
  plot(chl_time,vgpm[[i]],type='l',col='black',xlim=timedims,xaxt='n',yaxt='n'); 
  abline(v=ticks,lty=2)
  axis(1,labels=NA,at=ticks)
  axis(side=4)
  par(new=TRUE)	
  plot(chl_time,chl[[i]],type='l',col='dark green',xaxt='n',xlim=timedims,yaxt='n')
  mtext(adj=0,labs[i],cex=0.6)
  axis(1,labels=NA,at=ticks)
  mtext(paste(round(cor(chl[[i]],vgpm[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)	
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
timedims <- c(2003,2021)
cdims    <- c(5,30) 

#pdf('d:/dropbox/working/gradients/tzcf/plots/c_chl_timeseries.pdf',height=7,width=6)
pdf('~/dropbox/working/gradients/tzcf/plots/c_chl_timeseries_03_26_2021.pdf',height=6,width=8)
par(mfrow=c(6,1),mar=c(1,2,1,2),oma=c(2,3,2,3),cex.axis=0.7)

plot(chl_time,chl[[6]],type='l',col='dark green',yaxt='n',xlim=timedims,xaxt='n',ylim=c(0,3)); 
axis(side=2,col='dark green'); 
axis(1,labels=NA,at=ticks)
par(new=TRUE)
plot(chl_time,c[[6]],type='l',col='black',xaxt='n',yaxt='n',xlim=timedims)
axis(4)
abline(v=ticks,lty=2)
mtext(adj=0,labs[6],cex=0.6)
mtext(paste(round(cor(chl[[6]],c[[6]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)
for(i in 5:1){
  plot(chl_time,chl[[i]],type='l',col='dark green',yaxt='n',xlim=timedims,xaxt='n'); 
  axis(side=2,col='dark green'); 
  axis(1,labels=NA,at=ticks)
  
  par(new=TRUE)
  plot(chl_time,c[[i]],type='l',col='black',xaxt='n',yaxt='n',xlim=timedims)
  axis(4,col='black')
  abline(v=ticks,lty=2)
  mtext(adj=0,labs[i],cex=0.6)
  mtext(paste(round(cor(chl[[i]],c[[i]],use='pairwise.complete.obs'),3)),adj=1,cex=0.6)
}
axis(1,at=ticks)
mtext(outer=TRUE, side=2, expression('Chlorophyll Concentration [mg/m'^3*']'))
mtext(outer=TRUE, side=4, expression('Carbon Concentration [mg/m'^3*']'),line=0.75)
dev.off()


###########################################################
## LAG PLOTS ##############################################
###########################################################

lag=times_c=times_chl <- matrix(NA,nrow=17,ncol=6)

days <- interval(ymd('1990-01-01'),date_decimal(chl_time))
days <- as.period(days,unit='days')
days <- as.numeric(substring(days,1,4))

ii <- chl_time >2003 & chl_time <=2004
y0 <- 2003
y1 <- 2004
for(j in 1:6){
  for(i in 0:16){
    yy0 <- y0 + i
    yy1 <- y1 + i
    ii <- chl_time > yy0 & chl_time <= yy1
    t_chl <- chl_time[ii][chl[[j]][ii]==max(chl[[j]][ii],na.rm=TRUE) & !is.na(chl[[j]][ii])]
    t_c   <- chl_time[ii][c[[j]][ii]  ==max(c[[j]][ii],na.rm=TRUE)   & !is.na(c[[j]][ii])]
    lag[i+1,j] <- (t_chl-t_c)
  }
}

labs <- c(expression('00-10'*degree*'N'),expression('10-20'*degree*'N'),expression('20-30'*degree*'N'),
          expression('30-40'*degree*'N'),expression('40-50'*degree*'N'),expression('50-60'*degree*'N'))
yrs <- 2003:2019

muus=sdds <- numeric(6)
pdf('~/dropbox/working/gradients/tzcf/plots/lag_chl_c_03_29_2021.pdf',height=4,width=6)
  par(mfrow=c(1,1),cex.axis=0.8)
  plot(-999,xlim=c(0.5,6.5),ylim=c(-300,100),bty='n',xaxt='n',xlab='',ylab='')
  
  for(i in 1:6){
    muu <- mean(365*lag[,i])
    muus[i] <- muu
    sdd <- sd(365*lag[,i])
    sdds[i] <- sdd
    points(i,muu,pch=19,cex=1.2)
    segments(x0=i,x1=i,y0=muu-sdd,y1=muu+sdd)
  }
  abline(h=0,lty=2)
  axis(side=1,labels=labs,at=1:6)
  mtext(side=2,line=2.5,'Lag of Chlorophyll Maximum [Days]')
dev.off()


###########################################################
## CHL:C vs. MLD-averaged IRRADIANCE ######################
###########################################################
cols <-viridis(12)

pdf('~/dropbox/working/gradients/tzcf/plots/CHLC_vs_E.pdf',height=6,width=5)
par(mfrow=c(3,2),mar=c(2,1,1,1),oma=c(3,3,2,2),xpd=FALSE)
for(i in 6:1){ 
	plot(e[[i]],chlc[[i]],pch=19,cex=0.5,col=cols[month(ocnprd_date)],xlim=c(0,300),ylim=c(0,0.05))
	if(i==5){image.plot(matrix(1:12),col=cols,legend.only=TRUE)}
	mtext(labs[i],adj=0,cex=0.5)
}
mtext(side=1,outer=TRUE,expression('Mixed Layer Averaged Irradiance ['*mu*'E/m'^2*'/s]'),line=0.8) 
mtext(side=2,outer=TRUE,expression('Chlorophyll:Carbon Ratio [mgChl/mgC]'),line=1) 
dev.off()

##############################################################
## SCATTERPLOTES WITH MLD ####################################
##############################################################
mldin <- seq(0,150,0.1)

linep <- function(x,y,xin){
	fit <- lm(y ~ x)
	lines(mldin,predict(fit,newdata=list(x=xin)),col='dark grey')
}

# labs <- c('a) 50-60 Deg.N','g)','m)',
#           'b) 40-50 Deg.N','h)','n)',
#           'c) 30-40 Deg.N','i)','o)',
#           'd) 20-30 Deg.N','j)','p)',
#           'e) 10-20 Deg.N','k)','q)',
#           'f) 00-10 Deg.N','l)','r)')

labs <- c('a) 50-60 Deg.N','g)',
          'b) 40-50 Deg.N','h)',
          'c) 30-40 Deg.N','i)',
          'd) 20-30 Deg.N','j)',
          'e) 10-20 Deg.N','k)',
          'f) 00-10 Deg.N','l)')


pdf('~/dropbox/working/gradients/tzcf/plots/scatter_chl_c_npp_mld_03_26_2021.pdf',height=6,width=6)
par(mfrow=c(6,2),oma=c(3,2,2,2),mar=c(1,4,0,0),cex.axis=0.8,xpd=FALSE)
plot(mld[[6]],c[[6]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,70),xaxt='n',xlab='',ylab='')
linep(mld[[6]],c[[6]],mldin)
mtext(adj=0,cex=0.6,labs[1])
plot(mld[[6]],chl[[6]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,2),xaxt='n',xlab='',ylab='')
    linep(mld[[6]],chl[[6]],mldin)
    mtext(adj=0,cex=0.6,labs[2])
    # plot(mld[[6]],npp[[6]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,1200),xaxt='n',xlab='',ylab='')
	# 	linep(mld[[6]],npp[[6]],mldin)
	# 	mtext(adj=0,cex=0.6,labs[3])
		k = 3
for(i in 5:1){

	plot(mld[[i]],c[[i]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,40),xaxt='n',xlab='',ylab='')
		linep(mld[[i]],c[[i]],mldin)
		if(i==1) axis(1)
		if(i==3) mtext(side=2,expression('Carbon Concentration [mg/m'^3*']'),line=2,adj=0.335)
		mtext(adj=0,cex=0.6,labs[k])
		k=k+1
		plot(mld[[i]],chl[[i]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,0.5),xaxt='n',xlab='',ylab='')
		linep(mld[[i]],chl[[i]],mldin)
		if(i==1) axis(1)
		if(i==3) mtext(side=2,expression('Chlorophyll Concentration [mg/m'^3*']'),line=2,adj=0.335)
		mtext(adj=0,cex=0.6,labs[k])
		k=k+1
	# plot(mld[[i]],npp[[i]],pch=19,cex=0.25,xlim=c(0,150),ylim=c(0,1200),xaxt='n',xlab='',ylab='')
	# 	linep(mld[[i]],npp[[i]],mldin)
	# 	if(i==1) axis(1)
	# 	mtext(adj=0,cex=0.6,labs[k])
	# 	if(i==3) mtext(side=2,expression('Net Primary Productivity [mg/m'^2*'/day]'),line=2.0)
	# 	k=k+1
}
mtext(side=1,outer=TRUE,'Mixed Layer Depth [m]',line=1.5)
dev.off()
