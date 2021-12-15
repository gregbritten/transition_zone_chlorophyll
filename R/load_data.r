load(file='processed_data/PAR.rds') 
  PAR <- PAR*(1E6/86400)	#einstein/m2/day -> ueinstein/m2/s
load(file='processed_data/MLD.rds')	
load(file='processed_data/KD.rds')	
load(file='processed_data/I.rds')	
  I <- I*(1E6/86400)  #einstein/m2/day -> ueinstein/m2/s
load(file='processed_data/CHL.rds')	
load(file='processed_data/CHL_GIOP.rds')	
load(file='processed_data/CHL_GSM.rds')	
load(file='processed_data/BBP_GIOP.rds')	
load(file='processed_data/BBP_GSM.rds')	
load(file='processed_data/CARBON.rds')	
load(file='processed_data/VGPM.rds')	
load(file='processed_data/CBPM.rds')	
load(file='processed_data/EP_VGPM.rds')	
load(file='processed_data/CAFE.rds')	
load(file='processed_data/date.rds')	
load(file='processed_data/time.rds')	

load(file='processed_data/WOA_N.rds')
load(file='processed_data/WOA_depth.rds')

load(file='processed_data/ML_soda.rds')
load(file='processed_data/time_soda.rds')
soda_i <- 1:length(time_soda)
I_soda <- (1/(KD[,,soda_i]*ML_soda))*PAR[,,soda_i]*exp(-(KD[,,soda_i]*ML_soda))*(1-exp(-(KD[,,soda_i]*ML_soda)))


CARBON_GIOP <- 13000*(BBP_GIOP - 0.00035)
CARBON_GSM  <- 13000*(BBP_GSM  - 0.00035)

CHL_GSM[CHL_GSM>50] <- NA
CHL_GIOP[CHL_GIOP>50] <- NA

lat  <- seq(0,65,length.out=130)
lon  <- seq(-180,-115,length.out=130)
loni <- which(lon>=-180 & lon <=-115)

LAT <- list()
LAT[['lat00_10']] <- which(lat >=0  & lat <10)
LAT[['lat10_20']] <- which(lat >=10 & lat <20)
LAT[['lat20_30']] <- which(lat >=20 & lat <30)
LAT[['lat30_40']] <- which(lat >=30 & lat <40)
LAT[['lat40_50']] <- which(lat >=40 & lat <50)
LAT[['lat50_60']] <- which(lat >=50 & lat <60)
