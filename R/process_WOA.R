
library(ncdf4)

source('r/functions.r')
path <- "~/dropbox/data/woa/woa18/monthly/"

Nfiles <- list.files(paste0(path,"nitrate"))

nc    <- nc_open(paste0(path,"nitrate/",Nfiles[1]))
	lats  <- ncvar_get(nc,"lat")
	lons  <- ncvar_get(nc,"lon")
	woa_z <- ncvar_get(nc,"depth")
nc_close(nc)

lati <- which(lats>=0 & lats <= 65)
loni <- which(lons>=-180 & lons <=-115)

N <- array(NA,dim=c(130,130,43,12))

for(i in 1:12){
	ncN   <- nc_open(paste0(path,"nitrate/",Nfiles[i]))
	Ntmp  <- ncvar_get(ncN,"n_an")[loni,lati,]
	for(j in 1:43){
		N[,,j,i] <- resize_bilinear(z=Ntmp[,,j], yin=65, xin=65, yout=130, xout=130)	 
	}
}

save(file="processed_data/WOA_N.rds", N)
save(file="processed_data/WOA_depth.rds",woa_z)
