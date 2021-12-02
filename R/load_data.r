
#setwd('d:/dropbox/working/gradients/tzcf/data/')
setwd('~/dropbox/working/gradients/tzcf/data/')

library(lubridate)

lim <- 1:806

load('OCNPRD_CHL_gapfilled.rdata')
	CHL <- DAT[,,lim]
load('OCNPRD_MLD.rdata')
	MLD <- DAT
load('OCNPRD_carbon.rdata')
	C   <- DAT[,,lim]
load('OCNPRD_growth_rate.rdata')
	MU  <- DAT[,,lim]
load('OCNPRD_K490.rdata')
	K  <- DAT[,,lim]
load('OCNPRD_PAR.rdata')
	PAR<- (DAT[,,lim]/86400)*1E6  #einstein/m2/day -> ueinstein/m2/s 
load('OCNPRD_CAFE.rdata')
	CAFE<- CAFE[,,lim]
load('OCNPRD_VGPM.rdata')
	VGPM <- VGPM[,,lim]
load('OCNPRD_EPPLEY_VGPM.rdata')
	EPPLEY_VGPM <- EPPLEY_VGPM[,,lim]
load('OCNPRD_CBPM.rdata')
	CBPM <- CBPM[,,lim]
load('OCNPRD_BBP.rdata')
	BBP <- BBP[,,lim]

##############################################
## ALIGN OCNPRD AND MLD DATES ################
##############################################
load('OCNPRD_time.rdata')
load('OCNPRD_MLD_time.rdata')
mld_date    <- mld_time
ocnprd_date <- ocnprd_date[lim]
chl_time <- decimal_date(ymd(ocnprd_date))
mld_time <- decimal_date(ymd(mld_date))

MLD <- MLD[,,mld_date %in% ocnprd_date] #align the data arrays

Cint <- C*MLD



