require(tidyverse)

LLAtoECEF <- function(lat,long,alt){
  #conversion of geodetic lla to ECEF
  a <- 6378137 #equatorial radius in WGS-84
  f <- 1 / 298.257223563 #flatteningin WGS-84
  lat <- lat/180*pi
  long <- long/180*pi
  
  e <- (2*f-f^2)^0.5
  r <- a / (1 - e^2 * sin(lat) ^ 2)^0.5
 
  x <- (r + alt) * cos(lat) * cos(long)
  y <- (r + alt) * cos(lat) * sin(long)
  z <- (r * (1-e^2) + alt) * sin(lat)
  return(c(x=x,y=y,z=z))
}

ECEFtoLLA <- function(x,y,z){
  TOLERANCE <- 1e-6 #reuqired precision of longitude calculation
  a <- 6378137 #equatorial radius in WGS-84
  f <- 1 / 298.257223563 #flatteningin WGS-84
  e <- (2*f-f^2)^0.5
  b <- a*(1-f) 
  
  long <- atan2(y,x)
  p= (x^2+y^2)^0.5

  lat_ <- atan(z/p/(1-e^2))
  
  nextLat <- function(lat){
    r <- a / (1-e^2 *sin(lat)^2)^0.5
    h <- p/cos(lat)-r
    answer <- atan(z/p / (1- e^2 *r/(r+h) ))
    return(answer)
  }
  
  lat <- nextLat(lat_)
  while(abs(lat -lat_)>TOLERANCE){
    lat_ <- lat
    lat <- nextLat(lat)
  }
  
  h=p/cos(lat) - a / (1-e^2 *sin(lat)^2)^0.5

  long <- long*180/pi 
  lat <- lat*180/pi
   
  return(c(lat=lat,long=long,alt=h))
}

#This can be done in proj! EPSG:1165 is itrs 2014

ECEFtoENU <- function(lat,long){
  lat <- lat*pi/180
  long <- long*pi/180
  transformMatrix <- matrix(c(-sin(long),         cos(long)          ,0       ,
                              -cos(long)*sin(lat),-sin(long)*sin(lat),cos(lat),
                              cos(long)*cos(lat) , sin(long)*cos(lat),sin(lat)
                              ),byrow=TRUE,ncol=3,nrow=3)
  
  function(x,y,z){
    delta <- matrix(c(x,y,z),ncol=1) 
    print(delta)
    return(transformMatrix %*% delta)    
  }
}

getAzimuthELevation <- function(lat1,long1,alt1,lat2,long2,alt2){
  #reports wrt to device 
  pDevice <- LLAtoECEF(lat1,long1,alt1)
  pSat <- LLAtoECEF(lat2,long2,alt2)
  print(pSat)
  print(pDevice)
  delta <- pSat-pDevice
  print(delta)
  deltaENU <- ECEFtoENU(lat1,long1)(delta[["x"]],delta[["y"]],delta[["z"]]) 
  distance <- sum(delta^2)^0.5
  elevation <- asin(deltaENU[3,1]/distance)
  azimuth <- atan(deltaENU[1,1]/deltaENU[2,1])
  
  return(c(elevation=elevation,azimuth=azimuth))
}


