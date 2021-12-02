rm(list=ls())

library(rgdal)
library(fields)
library(viridis)
library(lubridate)

source('r/functions.r')

Nlon <- 720
Nlat <- 360

Nlat0 <- 1080
Nlon0 <- 2160

lons <- seq(-180,180,length.out=Nlon)
lats <- seq(-90,90,  length.out=Nlat)

loni <- which(lons >= -180 & lons <= -115)
lati <- which(lats >=    0 & lats <=   65)

nlon <- length(loni)
nlat <- length(lati)

lats_sv <- lats[lati]
lons_sv <- lons[loni]

save(file='processed_data/lat.rds',lats_sv)	
save(file='processed_data/lon.rds',lons_sv)	

#####################################################################################
## PROCESS VARIABLES ################################################################
#####################################################################################
raw_dir <- "/Users/gregorybritten/dropbox/data/ocean_productivity/"   #directory where raw HDF files are stored

par_files      <- list.files(paste0(raw_dir,'par'), full.names=TRUE)
mld_files      <- list.files(paste0(raw_dir,'mld'), full.names=TRUE)
kd_files       <- list.files(paste0(raw_dir,'k490'), full.names=TRUE)
chl_files      <- list.files(paste0(raw_dir,'modis_chl'), full.names=TRUE)
chl_giop_files <- list.files(paste0(raw_dir,'modis_chl_giop'), full.names=TRUE)
chl_gsm_files  <- list.files(paste0(raw_dir,'modis_chl_gsm'), full.names=TRUE)
bbp_giop_files <- list.files(paste0(raw_dir,'bbp_giop'), full.names=TRUE)
bbp_gsm_files  <- list.files(paste0(raw_dir,'bbp_gsm'), full.names=TRUE)
carbon_files   <- list.files(paste0(raw_dir,'cbpm_carbon'), full.names=TRUE)

nums <- substr(chl_giop_files,start = 74,80)[1:830]   #extract the complete date numbers for 2002:2020; chl_giop_files is complete where some others aren't

##construct dates
yr    <- as.numeric(substr(nums,1,4))
doy   <- as.numeric(substr(nums,5,7))
date  <- as.character(as.Date(doy, origin = paste(yr,"01","01",sep='-')))
time  <- decimal_date(ymd(date))

Nfiles <- length(nums)

I = PAR = MLD = KD = CHL = CHL_GIOP = CHL_GSM = BBP_GIOP = BBP_GSM = CARBON <- array(NA,dim=c(nlon,nlat,Nfiles))

for(i in 1:Nfiles){
  print(i)
  ii <- nums[i]

  par_file      <- paste0(raw_dir,'par/par.',ii,'.hdf')
  mld_file      <- paste0(raw_dir,'mld/mld.',ii,'.hdf')
  kd_file       <- paste0(raw_dir,'k490/k490.',ii,'.hdf')
  chl_file      <- paste0(raw_dir,'modis_chl/chl.',ii,'.hdf')
  chl_giop_file <- paste0(raw_dir,'modis_chl_giop/chl.',ii,'.hdf')
  chl_gsm_file  <- paste0(raw_dir,'modis_chl_gsm/chl.',ii,'.hdf')
  bbp_giop_file <- paste0(raw_dir,'bbp_giop/bbp.',ii,'.hdf')
  bbp_gsm_file  <- paste0(raw_dir,'bbp_gsm/bbp.',ii,'.hdf')
  carbon_file   <- paste0(raw_dir,'cbpm_carbon/carbon.',ii,'.hdf')

  if(par_file%in%par_files)           {par_tmp <- as.matrix(readGDAL(par_file))[,Nlat0:1]}           else{par_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(mld_file%in%mld_files)           {mld_tmp <- as.matrix(readGDAL(mld_file))[,Nlat0:1]}           else{mld_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(kd_file%in%kd_files)             {kd_tmp <- as.matrix(readGDAL(kd_file))[,Nlat0:1]}             else{kd_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(chl_file%in%chl_files)           {chl_tmp <- as.matrix(readGDAL(chl_file))[,Nlat0:1]}           else{chl_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(chl_giop_file%in%chl_giop_files) {chl_giop_tmp <- as.matrix(readGDAL(chl_giop_file))[,Nlat0:1]} else{chl_giop_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(chl_gsm_file%in%chl_gsm_files)   {chl_gsm_tmp <- as.matrix(readGDAL(chl_gsm_file))[,Nlat0:1]}   else{chl_gsm_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(bbp_giop_file%in%bbp_giop_files) {bbp_giop_tmp <- as.matrix(readGDAL(bbp_giop_file))[,Nlat0:1]} else{bbp_giop_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(bbp_gsm_file%in%bbp_gsm_files)   {bbp_gsm_tmp <- as.matrix(readGDAL(bbp_gsm_file))[,Nlat0:1]}   else{bbp_gsm_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  if(carbon_file%in%carbon_files)     {carbon_tmp <- as.matrix(readGDAL(carbon_file))[,Nlat0:1]}     else{carbon_tmp <- matrix(NA,nrow=Nlon0, ncol=Nlat0)}
  
  par_tmp[par_tmp==-9999]           <- NA
  mld_tmp[mld_tmp==-9999]           <- NA
  kd_tmp[kd_tmp==-9999]             <- NA
  chl_tmp[chl_tmp==-9999]           <- NA
  chl_giop_tmp[chl_giop_tmp==-9999] <- NA
  chl_gsm_tmp[chl_gsm_tmp==-9999]   <- NA
  bbp_giop_tmp[bbp_giop_tmp==-9999] <- NA
  bbp_gsm_tmp[bbp_gsm_tmp==-9999]   <- NA
  carbon_tmp[carbon_tmp==-9999]     <- NA
  
  kd_tmp[is.na(mld_tmp)]  <- NA
  par_tmp[is.na(mld_tmp)] <- NA
   
  PAR[,,i]      <- resize_bilinear(par_tmp,      xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  MLD[,,i]      <- resize_bilinear(mld_tmp,      xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  KD[,,i]       <- resize_bilinear(kd_tmp,       xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  CHL[,,i]      <- resize_bilinear(chl_tmp,      xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  CHL_GIOP[,,i] <- resize_bilinear(chl_giop_tmp, xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  CHL_GSM[,,i]  <- resize_bilinear(chl_gsm_tmp,  xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  BBP_GIOP[,,i] <- resize_bilinear(bbp_giop_tmp, xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  BBP_GSM[,,i]  <- resize_bilinear(bbp_gsm_tmp,  xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  CARBON[,,i]   <- resize_bilinear(carbon_tmp,   xin=Nlon0, xout=Nlon, yin=Nlat0, yout=Nlat)[loni,lati]
  I[,,i]        <- (1/(KD[,,1]*MLD[,,i]))*PAR[,,i]*exp(-(KD[,,i]*MLD[,,i]))*(1-exp(-(KD[,,i]*MLD[,,i])))
}

save(file='processed_data/PAR.rds',      PAR)	
save(file='processed_data/MLD.rds',      MLD)	
save(file='processed_data/KD.rds',       KD)	
save(file='processed_data/I.rds',        I)	
save(file='processed_data/CHL.rds',      CHL)	
save(file='processed_data/CHL_GIOP.rds', CHL_GIOP)	
save(file='processed_data/CHL_GSM.rds',  CHL_GSM)	
save(file='processed_data/BBP_GIOP.rds', BBP_GIOP)	
save(file='processed_data/BBP_GSM.rds',  BBP_GSM)	
save(file='processed_data/CARBON.rds',   CARBON)	
save(file='processed_data/date.rds',     date)	
save(file='processed_data/time.rds',     time)	


