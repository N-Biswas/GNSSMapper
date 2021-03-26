require(tidyverse)
require(polynom)

SP3_LOCAL <- "processedData/ephemeris"
IONEX_LOCAL <- "processedData/ionosphere"

#SP3 format is given here ftp://igs.org/pub/data/format/sp3c.txt
#ESA chosen data source as it records all constellations, unlike IGS
# metadata: http://navigation-office.esa.int/products/gnss-products/esm.acn
# NB: time is GPS time , ITRF2014 reference frame
# if the finals are not released, IGS rapids or ultras are used, but only gps/glonass available.

#RINEX MGXX basis http://mgex.igs.org/IGS_MGEX_Products.php
# https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/gnss_mgex.html

dateToGPS <- function(input){
  DATE_GPS_WEEK0 <- as.Date("1980-01-06")
  if (class(input)=="gpsDate") {
    return (input)
  } else {
    date<- as.Date(input)
    daysDiff <- as.numeric(date - DATE_GPS_WEEK0)
    gpsDate <- structure(c(week = daysDiff %/% 7,
                          day=daysDiff %% 7), 
                        class = "gpsDate")
    return(gpsDate)
  }
}

dateToYearDay <- function(input){
  date<- as.Date(input)
  year <- format(date,"%Y")
  YearDay <- c(year = format(date,"%y"),
               day = format(date,"%j"))
  return(YearDay)
}

getSP3File <- function(date){
  SP3_DATASITE <- "http://navigation-office.esa.int/products/gnss-products/"
  locateFile <- function(date){
    tryLocations <- function(local,url){
      if (!file.exists(local)){
        download.file(url,local)
      }
      return(local)
    }
    
    gpsDate <- dateToGPS(date)  
    filenameFinal <- paste0(c("esm",gpsDate["week"],gpsDate["day"],".sp3.gz"),collapse="")
    urlFinal <- paste0(c(SP3_DATASITE,gpsDate["week"],"/",filenameFinal),collapse="")
    localFinal <- paste0(c(getwd(),"/",SP3_LOCAL,"/",filenameFinal),collapse="")
    
    filenameRapid <- paste0(c("esr",gpsDate["week"],gpsDate["day"],".sp3.Z"),collapse="")
    urlRapid <- paste0(c(SP3_DATASITE,gpsDate["week"],"/",filenameRapid),collapse="")
    localRapid <- paste0(c(getwd(),"/",SP3_LOCAL,"/",filenameRapid),collapse="")
    
    filenameUltra <- paste0(c("esu",gpsDate["week"],gpsDate["day"],"_00.sp3.Z"),collapse="")
    urlUltra <- paste0(c(SP3_DATASITE,gpsDate["week"],"/",filenameUltra),collapse="")
    localUltra <- paste0(c(getwd(),"/",SP3_LOCAL,"/",filenameUltra),collapse="")
    
    if (!dir.exists(SP3_LOCAL)) {dir.create(SP3_LOCAL)}
    
    try(return(tryLocations(localFinal,urlFinal)))
    try(return(tryLocations(localRapid,urlRapid)))
    return(tryLocations(localUltra,urlUltra))
    }

  filename <- locateFile(date)
  message("Using ",filename)
  
  if (endsWith(filename,".gz")) {
    conn <- gzfile(filename)
    sp3 <- readLines(conn)
    close(conn)
  } else {
    sp3 <- system2("uncompress", args= c("-c",filename),stdout=TRUE)
  }

  return(sp3)
}

getOrbits <- function(sp3){
  #from an sp3 text file with format as given here ftp://igs.org/pub/data/format/sp3c.txt
  #returns tibble of locations for each epoch and svid
  KM_TO_METRES <- 1000
  MICRO_TO_NANO <- 1000
  
  getTime <- function(row){
    year <- substr(row,4,7)
    month <- substr(row,9,10)
    day <- substr(row,12,13)
    hour <- substr(row,15,16)
    min <- substr(row,18,19)
    sec <- substr(row,21,31)
    time <- c(year,"-",month,"-",day," ",hour,":",min,":",sec) %>% paste0(collapse="") %>% as.POSIXct(tz="UTC")
    return(time)
  }
  
  getXYZ <- function(row){     
    svid <- substr(row,2,4)
    x <- as.double(substr(row,5,18))*KM_TO_METRES
    y <- as.double(substr(row,19,32))*KM_TO_METRES
    z <- as.double(substr(row,33,46))*KM_TO_METRES
    clockError <- as.double(substr(row,47,60))*MICRO_TO_NANO
  return(list(svid=svid,x=x,y=y,z=z,clockError=clockError))}
  
  epochsTotal <- substr(sp3[[1]],33,39) %>% as.integer()  #grab number of epochs from the top line
  svidsTotal <- substr(sp3[[3]],4,6) %>% as.integer()   #count number of SVIDs
  out <- vector("list",epochsTotal*svidsTotal)
  time <- NULL
  i <- 1
  for (line in sp3) {
    if(startsWith(line,("* "))) {
      time <- getTime(line)
    } else if (startsWith(line,"P")) {
      out[[i]] <- c(list(time=time),getXYZ(line))
      i <- i+1
    }
    next()
  }
  
  out <- transpose(out) %>% simplify_all() %>%  as_tibble()
  class(out[["time"]])<- "POSIXct"
  return(out)
  
}

getSVIDLocation_Deprecated <- function(date){
  #this returns a function to estimate satellite position, having loaded and processed the ephemeris info for a single date
  #function assumes constant velocity between measured points. 

  orbits <- getSP3File(date) %>% getOrbits()
  orbits[["time"]] <- as.numeric(orbits[["time"]])
  times <- unique(orbits[["time"]]) %>% sort() 
  orbitsList <- orbits %>% split(.[["time"]]) %>% map(~.x[order(.x[["svid"]]),] %>% select(-time) %>% column_to_rownames("svid") %>% as.matrix())
  if (!all(map_lgl(orbitsList,~identical(rownames(.x),rownames(orbitsList[[1]]))))) {stop("The SP3 file is missing some SVs")}
  function(time){
    time <- as.numeric(time)
    preIndex <- findInterval(time,times,all.inside = TRUE)
    priorClock <- times[preIndex]
    postClock <- times[preIndex+1]
    weightPost <- (time-priorClock)/(postClock-priorClock)
    prior <- orbitsList[[as.character(priorClock)]]  
    post <- orbitsList[[as.character(postClock)]] 
    resultMatrix <- prior*(1-weightPost) + post*weightPost
    return(as_tibble(resultMatrix,rownames="Svid"))
  }
}

getSVIDLocation <- function(date){
  #given a uncorrected transmitter time and svid, this returns the corrected clock and satellite. 
  #the initial function generation loads and processes the ephemeris info for a single date.
  # Location is based on sp3 data and  uses a lagrange interpolation of a 7 order poly. This has an  s.d error of 3.3cm.following http://www.acc.igs.org/orbits/orbit-interp_gpssoln03.pdf (see also http://www.acc.igs.org/orbits/orbit-interp_gpssoln06.pdf)
  # clock error is based on sp3
  GPS_UTC_OFFSET <- as.integer(as.POSIXct('1980-01-06',tz="UTC"))
  NANOS_TO_SECOND <- 1e-9
  poly <-  function(i,allData){
    #solves a lagrange interpolation centred at index i. Values scaled to avoid calculation problems with small coefficients
    if(i<4 | i> nrow(allData) -4)(warning("outside fit interval"))
    data <- allData[(i-3):(i+4),]
    mid <- data[4,]
    scale <- (data[8,]-data[1,]) %>% replace(.==0,1)
    scaled_data <- sweep(data,2,mid) %>% sweep(2,scale,"/")
    polys_svid <- poly.calc(scaled_data[,1],scaled_data[,-1])
    return(list(scale=scale,mid=mid,polys_svid=polys_svid))
  }
  
  orbits <- getSP3File(date) %>% getOrbits()
  
  orbits[["time"]] <- (orbits[["time"]] %>% as.numeric()-GPS_UTC_OFFSET)
  orbits <- orbits[order(orbits[["time"]],orbits[["svid"]]),]
  svids <- unique(orbits[["svid"]]) #used to provide rownames for final dataframe
  orbitTimes <- unique(orbits[["time"]])
  l <- length(orbitTimes)
  orbitsXYZ <- orbits %>% split(.[["svid"]]) %>% map(~.x %>% select(-clockError,-svid) %>% as.matrix())# pivot_wider(id_cols=c(time,svid),values_from=c(x,y,z))
  orbitsT <- orbits %>% select(time,svid,clockError)%>% pivot_wider(id_cols=time,names_from=svid,values_from=clockError) %>% as.matrix()
  i <- seq(4, l -4, by=4)
  polysXYZ <- map(orbitsXYZ,~map(i,poly,.x))
  velocitiesXYZ <- map_depth(polysXYZ,2,~.x[["polys_svid"]] %>% deriv) 
  
  locateSatellite <- function(time,svid){
    if (!is.na(time) && svid %in%svids) {
      i <- findInterval(time,orbitTimes,all.inside = TRUE)
      if(i<3 | i> l-2)(warning("outside fit interval"))
      j <- (i+1) %/% 4
      poly <- polysXYZ[[svid]][[j]]
      predictVector <- (map_dbl(poly[["polys_svid"]],~predict(.x,(time-poly[["mid"]][1])/poly[["scale"]][1])) *poly[["scale"]][-1]+poly[["mid"]][-1])
    return(predictVector)
    } else {
      return(c(x=NA_real_,y=NA_real_,z=NA_real_))
    }
  }
  
  velocitySatellite <- function(time,svid){
    if (!is.na(time) && svid %in%svids)  {
      i <- findInterval(time,orbitTimes,all.inside = TRUE)
      if(i<3 | i> l-2)(warning("outside fit interval"))
      j <- (i+1) %/% 4
      poly <- polysXYZ[[svid]][[j]]
      velocity <- velocitiesXYZ[[svid]][[j]]
      predictVector <- (map_dbl(velocity,~predict(.x,(time-poly[["mid"]][1])/poly[["scale"]][1])) *poly[["scale"]][-1]/poly[["scale"]][1]) #using chain rule
      return(predictVector)
    } else {
      return(c(x=NA_real_,y=NA_real_,z=NA_real_))
    }
  }
  
  rotateECEF <- function(location,time){
    if(!is.na(time) && !is.na(location)){
      ROTATION_RATE <-7292115.0 * 10^-11 #WGS84 radians/second
      theta <- ROTATION_RATE*time
      transformMatrix <- matrix(c( cos(theta), sin(theta),0,
                                  -sin(theta), cos(theta),0,
                                   0         , 0         ,1)
                                ,byrow=TRUE,ncol=3,nrow=3)
        
      return(transformMatrix %*% location %>% drop() %>% setNames(c("x","y","z")))
    }
    else {
      return(c(x=NA_real_,y=NA_real_,z=NA_real_))
    }
  }
  
  getEphemerisClockError <- function(time,svid){
    #this is currently picking from prior 15 minute interval. I could get the proper 30second sample and use that.
    if (!is.na(time) && svid %in%svids) {
      return(orbitsT[findInterval(time*NANOS_TO_SECOND,orbitsT[,"time"],all.inside = TRUE),svid] %>% unname())
    } else {
      return(NA_real_)
    }
  }
  
  getRelativisticClockError <- function(location,velocity,svid){
    if(substr(svid,1,1)=="R"){
      return(0)
    } else {
      LIGHTSPEED <- 299792458
      dotproduct <- sum(location*velocity)
      return(-2/(LIGHTSPEED^2) * dotproduct/NANOS_TO_SECOND)
    }
  }
    
  estimateMeasurementTime <- function(receiverGPSTimeNanos,svid){
    #this estimates transmission time for a given epoch signal based on estimated satellite distance of 22,000 km  except for Beidou GEO/IGSO which are estimated at 38,000km
    # GPS_UTC_OFFSET <- as.integer(as.POSIXct('1980-01-06',tz="UTC"))
    # GPS_UTC_LEAPSECONDS <- -18
    LIGHTSPEED_NANOS <- 0.299792458
    BEIDOU_HIGH_SVID <- c("C01","C02","C03","C13","C16","C59","C31","C04","C05","C06","C07","C08","C09","C10","C38","C18","C39","C40")
    BEIDOU_HIGH <- 38000000/LIGHTSPEED_NANOS
    ORBIT <- 22000000/LIGHTSPEED_NANOS
    time <- receiverGPSTimeNanos #(as.numeric(EpochUTCTime)-GPS_UTC_OFFSET-GPS_UTC_LEAPSECONDS)*1e9
    return (time - ifelse (svid %in% BEIDOU_HIGH_SVID,BEIDOU_HIGH,ORBIT))
  }

  function(isEstimate,transmitterGPSTimeNanos,receiverGPSTimeNanos,svid){
    if (isEstimate){
      transmitterGPSTimeNanos <- estimateMeasurementTime(receiverGPSTimeNanos,svid) #approximate based on an average pseudorange
      locationECEFTransmitted <- locateSatellite(transmitterGPSTimeNanos*NANOS_TO_SECOND,svid) 
      return(c(time=transmitterGPSTimeNanos,locationECEFTransmitted,clockError=NA_real_,setNames(c(NA_real_,NA_real_,NA_real_),c("Rx","Ry","Rz")),setNames(c(NA_real_,NA_real_,NA_real_),c("Vx","Vy","Vz")),relativisticClockError=NA_real_))
    } else{
      clockError <- getEphemerisClockError(transmitterGPSTimeNanos,svid)
      transmitterGPSTimeNanos <- transmitterGPSTimeNanos-clockError  #first correct for clock offsets
      locationECEFTransmitted <- locateSatellite(transmitterGPSTimeNanos*NANOS_TO_SECOND,svid) #then locate the satellite
      locationECEFReceived <- rotateECEF(locationECEFTransmitted,(receiverGPSTimeNanos-transmitterGPSTimeNanos)*NANOS_TO_SECOND)
      rotation <- locationECEFReceived-locationECEFTransmitted
      velocityECEFTransmitted <- velocitySatellite(transmitterGPSTimeNanos*NANOS_TO_SECOND,svid) 
      RelativisticClockError <- getRelativisticClockError(locationECEFTransmitted,velocityECEFTransmitted,svid)
      transmitterGPSTimeNanos <- transmitterGPSTimeNanos-RelativisticClockError #correct for relativistic error
      return(c(time=transmitterGPSTimeNanos,locationECEFReceived,clockError=clockError,setNames(rotation,c("Rx","Ry","Rz")),setNames(velocityECEFTransmitted,c("Vx","Vy","Vz")),relativisticClockError=RelativisticClockError))
    }
  }
  
}

getIonexFile <- function(date){
  IONEX_DATASITE <- "http://navigation-office.esa.int/products/gnss-products/"
  locateFile <- function(date){
    tryLocations <- function(local,url){
      if (!file.exists(local)){
        download.file(url,local)
      }
      return(local)
    }
    
    gpsDate <- dateToGPS(date)  
    yearDayDate <- dateToYearDay(date)  
    filenameFinal <- paste0(c("esag",yearDayDate["day"],"0.",yearDayDate["year"],"i.Z"),collapse="")
    urlFinal <- paste0(c(IONEX_DATASITE,gpsDate["week"],"/",filenameFinal),collapse="")
    localFinal <- paste0(c(getwd(),"/",IONEX_LOCAL,"/",filenameFinal),collapse="")
    
    filenameRapid <- paste0(c("esrg",yearDayDate["day"],"0.",yearDayDate["year"],"i.Z"),collapse="")
    urlRapid <- paste0(c(IONEX_DATASITE,gpsDate["week"],"/",filenameRapid),collapse="")
    localRapid <- paste0(c(getwd(),"/",IONEX_LOCAL,"/",filenameRapid),collapse="")
    
    filenameUltra <- paste0(c("ehrg",yearDayDate["day"],"0.",yearDayDate["year"],"i.Z"),collapse="")
    urlUltra <- paste0(c(IONEX_DATASITE,gpsDate["week"],"/",filenameUltra),collapse="")
    localUltra <- paste0(c(getwd(),"/",IONEX_LOCAL,"/",filenameUltra),collapse="")
    
    if (!dir.exists(IONEX_LOCAL)) {dir.create(IONEX_LOCAL)}
    
    try(return(tryLocations(localFinal,urlFinal)))
    try(return(tryLocations(localRapid,urlRapid)))
    return(tryLocations(localUltra,urlUltra))
  }
  
  filename <- locateFile(date)
  message("Using ",filename)
  
  if (endsWith(filename,".gz")) {
    conn <- gzfile(filename)
    ionex <- readLines(conn)
    close(conn)
  } else {
    ionex <- system2("uncompress", args= c("-c",filename),stdout=TRUE)
  }
  
  return(ionex)
}

getHeights <- function(ionex) {
  KILOMETRES_TO_METRES=1000
  iono <- substr(ionex[[16]],3,8) %>% as.integer() * KILOMETRES_TO_METRES
  earth <- substr(ionex[[14]],3,8) %>% as.integer() * KILOMETRES_TO_METRES
  return(c(iono=iono,earth=earth))
}

getTECgrid <- function(ionex){
  #from an ionex text file with format as given here ftp://www.igs.org/pub/data/format/ionex1.pdf
  #assumes 2d map.
  #returns tibble of TEC for each epoch and latlonggrid
  
  getTime <- function(row){
    year <- substr(row,1,6)
    month <- substr(row,7,12)
    day <- substr(row,13,18)
    hour <- substr(row,19,24)
    min <- substr(row,25,30)
    sec <- substr(row,31,36)
    time <- c(year,"-",month,"-",day," ",hour,":",min,":",sec) %>% paste0(collapse="") %>% as.POSIXct()
    return(time)
  }
  
  getLat <- function(row){
    return(substr(row,3,8) %>% as.numeric())
  }
  getTEC <- function(index){     
    start <- index+1
    end <- index+ceiling(length(longRange)/16)
    readings <- paste0(ionex[start:end],collapse="")
    chop <- seq(1,length(longRange)*5,by=5)
    TEC <- map(chop,~(substr(readings,.x,.x+5) %>% as.double())*exponent)
    names(TEC) <- longRange
    return(TEC)
  }
  
  exponent <- 10^(substr(ionex[[19]],1,6) %>% as.integer())
  epochsTotal <- substr(ionex[[8]],1,6) %>% as.integer()  #grab number of epochs
  lat1 <- substr(ionex[[17]],3,8) %>% as.numeric() 
  lat2 <- substr(ionex[[17]],9,14) %>% as.numeric()
  dlat  <- substr(ionex[[17]],15,20) %>% as.numeric()
  latTotal <- ((lat2-lat1)/ dlat) %>% as.integer() +1 #count number of latitudes
  long1 <- substr(ionex[[18]],3,8) %>% as.numeric() 
  long2 <- substr(ionex[[18]],9,14) %>% as.numeric()
  dlong  <- substr(ionex[[18]],15,20) %>% as.numeric() 
  longRange <- seq(long1,long2,by=dlong)
  out <- vector("list",epochsTotal*latTotal)
  time <- NULL
  i <- 1
  start <- detect_index(ionex,endsWith,"START OF TEC MAP    ")
  end <- detect_index(ionex,endsWith,"START OF RMS MAP    ")
  for (n in seq(start,end,by=1)) {
    if(endsWith(ionex[[n]],"EPOCH OF CURRENT MAP")) {
      time <- getTime(ionex[[n]])
    } else if (endsWith(ionex[[n]],"LAT/LON1/LON2/DLON/H")) {
      out[[i]] <-c(list(time=time,lat=getLat(ionex[[n]])),getTEC(n))
      i <- i+1
    }
    next()
  }
 out <- transpose(out) %>% simplify_all() %>%  as_tibble()
 class(out[["time"]])<- "POSIXct"
 return(out)
}

getTEC <- function(grid){
  GPS_UTC_OFFSET <- as.integer(as.POSIXct('1980-01-06',tz="UTC"))
  NANOS_TO_SECOND <- 1e-9
  grid[["time"]] <- grid[["time"]] %>% as.numeric()-GPS_UTC_OFFSET*NANOS_TO_SECOND
  grid <- grid[order(grid[["time"]],grid[["lat"]]),]
  lats <- unique(grid[["lat"]]) #used to provide rownames for final dataframe
  longs <- names(grid[c(-1,-2)]) %>% as.numeric()
  times <- unique(grid[["time"]])
  TECgrid <- grid %>% split(.[["time"]]) %>%  map(~.x %>% select(-time) %>% column_to_rownames("lat") %>% as.matrix())
  rotate <- function(DeltaNanos){
    #there is a sstrong correlation betweeon ionosphere and sun position. We compensate for the time mismatch by rotating map 
    NANOS_IN_DAY <- 1e9*3600*24
    return(DeltaNanos/NANOS_IN_DAY*360)
  }
  
  bound <- function(x,vec){
    x <- max(x,min(vec))
    x <- min(x,max(vec))
    return(x)
  }
  interpolate <- function(matrix,lat_,long_){
    i <-  findInterval(lat_,lats,all.inside=TRUE)
    j <- findInterval(long_,longs,all.inside=TRUE)
    p <- (lat_-lats[i])/(lats[i+1]-lats[i])
    q <- (long_-longs[j])/(longs[j+1]-longs[j])
    return(   (1-p) * (1-q) * matrix[i,  j  ] +
              p     * (1-q) * matrix[i+1,j  ] +
              (1-p) * q     * matrix[i  ,j+1] +
                p   * q     * matrix[i+1,j+1]
    )
    
  }
  
  function(lat,long,GPSTimeNanos){
    lat <- bound(lat,lats)
    long <- bound(long,longs)
    GPSTimeNanos <- bound(GPSTimeNanos,times)
    i <- findInterval(GPSTimeNanos,times,all.inside=TRUE)
    p <- (GPSTimeNanos-times[i])/(times[i+1]-times[i])
    prior <- interpolate(TECgrid[[i]],lat,long + rotate(GPSTimeNanos-times[i]))
    post <-  interpolate(TECgrid[[i+1]],lat,long + rotate(GPSTimeNanos-times[i+1]))
    return(prior * (1-p) + post * p)
  }
}

getIonosphericDelay <-  function(date){
  ionex <- getIonexFile(date)
  TECgrid <-  getTECgrid(ionex)
  heights <-getHeights(ionex)
  interpolateTEC <- getTEC(TECgrid)
  remove(ionex)
  remove(TECgrid)
  freq <- 1575.42*1e6
  alpha <-40.3e16/(freq^2)
  
  ionoHeight <- heights["iono"]
  RADIUS_EARTH <- heights["earth"]
  function (elevation,azimuth,latitude,longitude,GPSTimeNanos) {
  
    
    
    
    
    
    
    
    latRad <- latitude / 360 * 2 * pi
    longRad <- longitude / 360 * 2 * pi
    elevationRad <- elevation / 360 * 2 * pi
    azimuthRad <- azimuth / 360 * 2 * pi
    
    theta <- pi/2 - elevationRad - asin(RADIUS_EARTH/(RADIUS_EARTH+ionoHeight) * cos(elevationRad))
    latIPPRad <- asin(sin(latRad) * cos(theta) + cos(latRad)* sin(theta) * cos(azimuthRad) )
    longIPPRad <- longRad + theta * sin(azimuthRad)/cos(latIPPRad)
    latIPP <- latIPPRad * 360/2/pi
    longIPP <- longIPPRad *360/2/pi
    VTEC <- interpolateTEC(latIPP,longIPP,GPSTimeNanos)
     
    slant <- (1-(RADIUS_EARTH/(RADIUS_EARTH+ionoHeight) * cos(elevationRad))^2)^-0.5
    return(alpha * slant * VTEC)
  }
}

getTroposphericDelay <- function(){
  K_1 <- 77.604 #K/mbar
  K_2 <- 382000 #K^2/mbar 
  R_D <- 287.054 #J/(kg K)
  G_M <- 9.784 #m/s^2
  G <- 9.80665 #m/s^2
  paramAvgTable <- tribble(
    ~lat, ~pressure, ~temp,  ~vapourPressure, ~tempLapse, ~vapourLapse,
    15,   1013.25,   299.65, 26.31,           6.30e-3,    2.77,
    30,   1017.25,   294.15, 21.79,           6.05e-3,    3.15,
    45,   1015.75,   283.15, 11.66,           5.58e-3,    2.57,
    60,   1011.75,   272.15, 6.78,            5.39e-3,    1.81,
    75,   1013.00,   263.65, 4.11,            4.53e-3,    1.55
  ) %>% as.matrix()
  
  paramDeltaTable <- tribble(
    ~lat, ~pressure, ~temp,  ~vapourPressure, ~tempLapse, ~vapourLapse,
    15,   0,         0,      0,               0,          0,
    30,   -3.75,     7,      8.85,            0.25e-3,    0.33,
    45,   -2.25,     11,     7.24,            0.32e-3,    0.46,
    60,   -1.75,     15,     5.36,            0.81e-3,    0.74,
    75,   -0.5,      14,     3.39,            0.62e-3,    0.30
  ) %>% as.matrix()
  
  interpolateTable <- function(table,latitude){
    absLat <- abs(latitude)
    i <- findInterval(absLat,table[,"lat"])
    if (i==0){
      t <- table[1,-1]
    } else if(i==5){
      t <- table[5,-1]
    } else {
      p <- (absLat-table[i,"lat"])/15
      t <- table[i,-1]*(1-p) + table[i+1,-1]*(p)
    } 
    return(t)
  }
  
  function(elevation,height,latitude,date){
  day <- dateToYearDay(date)["day"] %>% as.numeric()
  elevationRad <- elevation / 360 * 2 * pi
  
  paramAvg <- interpolateTable(paramAvgTable,latitude)
  paramDelta <- interpolateTable(paramDeltaTable,latitude)
  dayMin <- if (latitude>0) 28 else 211
  param <- paramAvg-paramDelta * cos(2*pi/365.25*(day-dayMin))
  
  dryDelay0 <- 1e-6 * K_1 * R_D / G_M * param["pressure"]
  wetDelay0 <- 1e-6 * K_2 * R_D / ((param["vapourLapse"] + 1) * G_M - (param["tempLapse"] * R_D)) * param["vapourPressure"] / param["temp"]
  dryDelay <- (1-(param["tempLapse"] * height / param["temp"])) ^ (G / R_D / param["tempLapse"]) * dryDelay0
  wetDelay <- (1-(param["tempLapse"] * height / param["temp"])) ^ (((param["vapourLapse"]+1) *G / R_D / param["tempLapse"])-1) * wetDelay0
  obliquity <- 1.001 / (0.002001 + sin(elevationRad)^2)^0.5
  return((dryDelay+wetDelay)*obliquity)
  }
}

