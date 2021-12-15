##--put a clean line on a scatterplot--#####################
linep <- function(x,y,xin){
  fit <- lm(y ~ x)
  lines(mldin,predict(fit,newdata=list(x=xin)),col='dark grey')
}

##--modified image plot
image.plot2 <- function(x,y,z,zlim,col){
  tmp <- z
  tmp[tmp<zlim[1]] <- zlim[1]
  tmp[tmp>zlim[2]] <- zlim[2]
  image.plot(x,y,tmp,zlim=zlim,col=col,xlab='',ylab=''); box()
}

image2 <- function(x,y,z,zlim,col){
  tmp <- z
  tmp[tmp>zlim[2]] <- zlim[2]
  tmp[tmp<zlim[1]] <- zlim[1]
  image(x,y,tmp,col=col)
}

##--compute annual mean climatology
mapmean <- function(X,mnth){
  tmp <- X[,,month(date)==mnth]
  return(apply(tmp,c(1,2),function(x) mean(x,na.rm=TRUE)))
}

##--plot climatology
plotmap <- function(x,y,z,zlim,lab1,lab2,lab3,legend=TRUE){
  image2(x=x,y=y,z=z,col=viridis(25),zlim=zlim)
  map(xlim=c(-180,-110),ylim=c(0,70),col='grey',fill=TRUE,add=TRUE) 
  box()
  mtext(adj=0,lab1,cex=0.8)
  mtext(adj=-0.1,lab2)
  mtext(adj=0,lab3,line=1.5)
  if(legend==TRUE){image.plot(matrix(zlim),legend.only=TRUE,col=viridis(25))}
}

##--bilinear interpolation
resize_bilinear <- function(xin,yin,xout,yout,z){
  library(fields)
  obj   <- list(x=1:xin, y=1:yin, z = z)
  tempx <- seq(1,xin,length.out=xout)
  tempy <- seq(1,yin,length.out=yout)
  loc   <- make.surface.grid(list(tempx,tempy))
  look  <- interp.surface(obj,loc)
  return(as.surface(loc,look)$z)
}

##--COMPUTE AVERAGE OF LAT/LON SQUARE--##########################
latmean <- function(X, loni, lati){
  tmp <- X[loni,lati,]
  return(apply(tmp,3,function(x) mean(x,na.rm=TRUE)))
}

##--COMPUTE AVERAGE OF RATIO FOR LAT/LON SQUARE--##########################
latmeanratio <- function(X,Y,loni,lati){
  tmpx <- X[loni,lati,]
  tmpy <- Y[loni,lati,]
  meanx <- apply(tmpx,3,function(x) mean(x,na.rm=TRUE))
  meany <- apply(tmpy,3,function(x) mean(x,na.rm=TRUE))
  return(meanx/meany)
}


##--LAG OF YEARLY MAXIMUM--###########################
findlag <- function(X,Y){
  lag <- matrix(NA,nrow=17,ncol=6)
  ii <- time >2003 & time <=2004
  y0 <- 2003
  y1 <- 2004
  for(j in 1:6){
    xtmp <- latmean(X,loni=loni,lati=LAT[[j]])
    ytmp <- latmean(Y,loni=loni,lati=LAT[[j]])
    for(i in 0:16){
      yy0 <- y0 + i
      yy1 <- y1 + i
      ii <- time > yy0 & time <= yy1
      t_chl <- time[ii][xtmp[ii]  ==max(xtmp[ii],na.rm=TRUE)   & !is.na(xtmp[ii])]    #find the time of year when chlorophyll is maximum
      t_c   <- time[ii][ytmp[ii]  ==max(ytmp[ii],na.rm=TRUE)   & !is.na(ytmp[ii])]
      lag[i+1,j] <- (t_chl-t_c)
    }
  }
  return(lag)
}

##--PLOT LAG OF YEARLY MAXIMUM--#######################
plotlags <- function(lag){
  means <- apply(lag,2,mean)*365
  sds   <- apply(lag,2,sd)*365
  points(1:6,means,pch=19,cex=1.2)
  segments(x0=1:6,x1=1:6,y0=means-sds,y1=means+sds)
}


##--COMPUTE MONTHLY MEANS--#################################
monthmean <- function(X){
  a <- matrix(NA,nrow=6,ncol=12)
  for(i in 1:6){
    for(j in 1:12){
      b      <- X[,LAT[[i]],mnth==j]
      a[i,j] <- mean(b,na.rm=TRUE)
    }
  }
  return(a)
}


##--CORRELATION MAP--#####################################
cormap <- function(X,Y){
  tmp_cor <- matrix(NA,nrow=130,ncol=130)
  for(i in 1:130){
    for(j in 1:130){
      if(sum(!is.na(X[i,j,]))>5 & sum(!is.na(Y[i,j,]))>5)
        tmp_cor[i,j] <- cor(X[i,j,],Y[i,j,],use='pairwise.complete.obs')
    }
  }
  return(tmp_cor)
}

##--PLOT CORRELATION MAP--################################
plotcormap <- function(Z,lab){
  image(x=lon,y=lat,Z,zlim=c(-0.8,0.8),axes=FALSE,col=brewer.pal(n=9,name='RdBu'),xlab='',ylab='')
  map(xlim=c(-180,-115),ylim=c(0,65),fill=TRUE,col='grey',add=TRUE)
  box()
  axis(side=1); axis(side=2); mtext(side=c(3),c(lab),adj=0)
  abline(v=seq(-180,-115,10),lty=3); abline(h=seq(0,60,10),lty=3)
}

##--COMPUTE ZONAL CORRELATIONS--############################
zonalcor <- function(cormap){
  return(apply(cormap,2,mean,na.rm = TRUE))
}

##--INTERPOLATE SODA FIELDS--################################
interp <- function(X,timex,timey){
  Y <- array(NA,dim=c(130,130,length(timey)))
  for(i in 1:130){
    for(j in 1:130){
      tmp <- X[i,j,]
      if(sum(!is.na(tmp))>10){
        Y[i,j,] <- approx(x=timex,y=tmp,xout=timey)$y
      }
    }
  }
  return(Y)
}

