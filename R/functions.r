
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



