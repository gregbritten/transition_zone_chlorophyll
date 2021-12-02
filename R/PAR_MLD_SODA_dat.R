
rm(list=ls())
library(viridis)
library(fields)

#source('d:/dropbox/working/gradients/tzcf/github/greg/load_data.r')
source("~/dropbox/working/gradients/tzcf/github/greg/load_data.r")

load('~/dropbox/working/gradients/tzcf/data/ML.rdata')
load('~/dropbox/working/gradients/tzcf/data/SODA_time.rdata')

I <- array(NA,dim=c(130,130,806))

for(i in 1:806){
print(i)  
  par <- PAR[,,i]
  k   <- K[,,i]
  mld <- MLD[,,i]
  I[,,i]   <-   (1/(k*mld))*par*exp(-(k*mld))*(1-exp(-(k*mld)))
}

save(file='~/dropbox/working/gradients/tzcf/data/PAR_MLD.rdata',I)


###################################################################

ocnprd_time <- decimal_date(ymd(ocnprd_date))  #convert ocnprd to decimal year
soda_clip   <- time>ocnprd_time[1] & time<ocnprd_time[length(ocnprd_time)]  #subset SODA for OCNPRD times
SODA_ML     <- ML[,,soda_clip]
soda_time   <- time[soda_clip]

ocnprd_clip <- ocnprd_time>soda_time[1] & ocnprd_time<soda_time[length(soda_time)] #subset OCNPRD for SODA times
ocnprd_time_clip <- ocnprd_time[ocnprd_clip]
ocnprd_date_clip <- ocnprd_date[ocnprd_clip]

MLD_clip <- MLD[,,ocnprd_clip]
PAR_clip <- PAR[,,ocnprd_clip]
K_clip   <- K[,,ocnprd_clip]

SODA_ML_clip <- array(NA,dim=c(130,130,711))

for(i in 1:130){
  for(j in 1:130){
    soda_ts <- SODA_ML[i,j,]
    if(sum(!is.na(soda_ts))>10){
    afun <- approxfun(x=soda_time,y=soda_ts)
    SODA_ML_clip[i,j,] <- afun(ocnprd_time_clip)
    }
  }
}


PAR_MLD_SODA <- array(NA,dim=c(130,130,711))

for(i in 1:711){
  print(i)  
  par <- PAR_clip[,,i]
  k   <- K_clip[,,i]
  #mld <- MLD_clip[,,i]
  mld <- SODA_ML_clip[,,i]
  PAR_MLD_SODA[,,i]   <-   (1/(k*mld))*par*exp(-(k*mld))*(1-exp(-(k*mld)))
}

save(file='~/dropbox/working/gradients/tzcf/data/PAR_MLD_SODA.rdata',PAR_MLD_SODA)
save(file='~/dropbox/working/gradients/tzcf/data/OCNPRD_time_clipped.rdata',ocnprd_date_clip)




