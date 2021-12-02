
library(viridis)
library(fields)
library(lubridate)

#load(file="d:/dropbox/working/gradients/tzcf/data/WOA_N.rdata")
#load(file="d:/dropbox/working/gradients/tzcf/data/WOA_P.rdata")
#load(file="d:/dropbox/working/gradients/tzcf/data/WOA_depth.rdata")
load(file="~/dropbox/working/gradients/tzcf/data/WOA_N.rdata")
load(file="~/dropbox/working/gradients/tzcf/data/WOA_P.rdata")
load(file="~/dropbox/working/gradients/tzcf/data/WOA_depth.rdata")

#load("d:/dropbox/working/gradients/tzcf/data/ocnprd_MLD.rdata")
load("~/dropbox/working/gradients/tzcf/data/ocnprd_MLD.rdata")
MLD <- DAT
#load("d:/dropbox/working/gradients/tzcf/data/ocnprd_PAR.rdata")
load("~/dropbox/working/gradients/tzcf/data/ocnprd_PAR.rdata")
PAR <- DAT
#load("d:/dropbox/working/gradients/tzcf/data/ocnprd_k490.rdata")
load("~/dropbox/working/gradients/tzcf/data/ocnprd_k490.rdata")
K490 <- DAT
load("~/dropbox/working/gradients/tzcf/data/OCNPRD_time.rdata")
load("~/dropbox/working/gradients/tzcf/data/OCNPRD_MLD_time.rdata")
time <- decimal_date(ymd(ocnprd_date))
mnth_mld <- month(mld_time)

mnth <- month(ocnprd_date)

js <- c(1,5,10,20,30,45,55,60,70,80,90,100)
i  <-1:110
zz <- 1:30
zs <- z[zz]
lats <- seq(0,65,length.out=130)

#################################################################
js <- c(1,5,10,20,30,45,55,60,70,80,90,95)
js <- c(50,55,60,65,70,75,80,85,90)
js <- round(seq(50,90,length.out=12))
	
ZMLD=DN=PARS=K=E <- matrix(NA,12,length(js))
latts <- numeric(length(js))

for(l in 1:length(js)){
	j <- js[l]
	latts[l] <- lats[j] 
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
}


#par(mfrow=c(2,1),mar=c(2,2,2,4))
#matplot(ZMLD,type='l',lty=1,col=cols[1:12])
#image.plot(legend.only=TRUE,as.matrix(latts),col=cols[1:12])
#matplot(DN,type='l',lty=1,col=cols[1:12])

DMLD <- apply(rbind(ZMLD[12,],ZMLD),2,diff)

#matplot(DMLD*DN,type='l',lty=1,col=cols[1:12])
#image.plot(legend.only=TRUE,as.matrix(latts),col=cols[1:12])

cols <- viridis(13)

is <- c(4:12,1:3)

#############################################################################
## PLOT EACH 
par(mfrow=c(5,3),mar=c(2,4,0,0),oma=c(2,2,2,2))
matplot(ZMLD[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='')
	mtext(side=2,'MLD',line=2.5)
matplot(DMLD[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='')
	mtext(side=2,'dMLD/dt',line=2.5)
	abline(h=0,lty=2)

matplot(DN[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='',ylim=c(-0.5,4))
abline(h=0,lty=2)
	mtext(side=2,expression(italic('N'[0]*' - N')),line=2.5)
for(i in 1:length(js)){
	plot(apply(as.matrix(DMLD[is,i]*DN[is,i]),1,function(x) max(c(0,x))),
		type='l',lty=1,col=cols[i],ylim=c(0,60),bty='n',xaxt='n',xlab='',ylab='')
	mtext(round(lats[js[i]]),cex=0.8,adj=0.1,line=-2)
	#abline(h=0,lty=2)
	if(i%in%c(1,4,7)) mtext(side=2,expression(italic("w'N'")),line=2.5)
	if(i%in%c(10,11,12)) axis(side=1,at=1:12,labels=is)
	mtext(side=1,outer=TRUE,'Month',line=0.5)
}

CairoPDF('d:/dropbox/working/gradients/tzcf/plots/entraintment_flux_5x1.pdf',height=8,width=7)
#pdf('d:/dropbox/working/gradients/tzcf/plots/entraintment_flux_5x1.pdf',height=8,width=7,useDingbats=FALSE)
par(mfrow=c(5,1),mar=c(1,4,0,8),oma=c(4,2,2,2),xpd=FALSE,cex.axis=0.65)
	matplot(ZMLD[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='')
		mtext(side=2,'MLD',line=3.5); 		mtext(side=2,'[m]',line=2.5,cex=0.7)

		image.plot(legend.only=TRUE,as.matrix(lats[js]),col=cols,legend.shrink=1,legend.lab='Latitude')
		matplot(DMLD[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='')
		mtext(side=2,'dMLD/dt',line=3.5)
		mtext(side=2,'[m / month]',line=2.25,cex=0.7)
		abline(h=0,lty=2)

	matplot(DN[is,],type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='',ylim=c(-0.5,4.5))
	abline(h=0,lty=2)
		mtext(side=2,expression(italic('N'[0]*' - N')),line=3.5)
		mtext(side=2,expression('[mmol N / m'^3*']'),line=2.25,cex=0.7)

	plot(-999,ylim=c(0,60),xlim=c(1,12),type='n',xaxt='n',yaxt='n',bty='n',ylab='',xlab='')
	for(i in 1:length(js)){
		lines(apply(as.matrix(DMLD[is,i]*DN[is,i]),1,function(x) max(c(0,x))),
			lty=1,col=cols[i])
		#mtext(round(lats[js[i]]),cex=0.8,adj=0.1,line=-2)
		#abline(h=0,lty=2)
		#if(i%in%c(1,4,7)) mtext(side=2,expression(italic("w'N'")),line=2.5)
		#mtext(side=1,outer=TRUE,'Month',line=0.5)
	}
		axis(side=2)
		mtext(side=2,expression(italic("w'N'")),line=3.5)
		mtext(side=2,expression('[mmol N / m'^2*' / month]'),line=2.25,cex=0.7)

	matplot(log10(E[is,]),type='l',lty=1,col=cols,bty='n',xaxt='n',xlab='',ylab='',ylim=c(-2.5,1.5))
		mtext(side=2,expression(italic(bar(PAR)['z'])),line=3.5)
		mtext(side=2,expression('[log10('*mu*'mol/m'^2*')]'),line=2.25,cex=0.7)
		axis(side=1,at=1:12,labels=is,cex.axis=0.9)
		
	mtext(side=1,outer=TRUE,'Month',line=1,cex.axis=1)
dev.off()



