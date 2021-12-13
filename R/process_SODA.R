#longitudes are given 0-360

rm(list=ls())
library(ncdf4)

load(file='processed_data/time.rds')	
files_all <- list.files('~/dropbox/data/soda/five_day/',full.names=TRUE)[1600:2774]

##--GET LATS AND LONS OF SODA--########################
nc <- nc_open(files_all[1])
	lats <- ncvar_get(nc,'yt_ocean')
	lons <- ncvar_get(nc,'xt_ocean') - 180
nc_close(nc)

##--GET TIME OF SODA--################################
soda_time_all <- numeric(length(files_all))
for(i in 1:length(files_all)){
  print(i)
  nc               <- nc_open(files_all[i])
  soda_time_all[i] <- ncvar_get(nc,'time') 
  nc_close(nc)
}
soda_time_all <- (soda_time_all - 7307)/365 + 2000  #convert to decimal year

##--SUBET SODA FILES FOR TIME WITH BUFFER--#####################
files    <- files_all[soda_time_all >= time[1]-0.25]  #subset the soda_times for those greater or equal to the earliest ocnprd time, with a half year buffer for the interpolation
soda_time <- soda_time_all[soda_time_all >= time[1]-0.25]

##--SUBSET OCNPRD TIMES TO INTERPOLATE ONTO--###################
time_soda <- time[time <= soda_time[length(soda_time)]]


##-PROCESS DATA--#############################################
loni         <- which(lons> -180 & lons< -115)
lati         <- which(lats>    0 & lats< 65)
lat          <- lats[lati]
lon          <- lons[loni]

lonshift     <- c(360:720,1:359)
ML <- array(NA,dim=c(length(loni),length(lati),length(files)))
for(i in 1:length(files)){
	print(i)
  nc       <- nc_open(files[i])
	tmpML    <- ncvar_get(nc,'mlp')[lonshift,]
	ML[,,i]  <- tmpML[loni,lati]
	nc_close(nc)
}

##--INTERPOLATE TO OCEAN PROD DATABASE TIME AXIS--######################
ML_soda <- interp(ML,timex=soda_time,timey=time_soda)

##--SAVE FILES--#####################################
save(file='processed_data/ML_soda.rds',ML_soda)
save(file='processed_data/time_soda.rds',time_soda)






