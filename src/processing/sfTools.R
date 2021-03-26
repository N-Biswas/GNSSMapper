#tools to turn logged GNSS data into sf format
require(sf)
require(tidyverse)
source("code/processing/sp3.R")

CRS_LONGLAT_WGS84 <- 4979
CRS_CART_WGS84 <- "+proj=cart +datum=WGS84 +no_defs +wktext"
CRS_BNG <- 27700
ORIGIN <- c(x=3980000 ,y=-10000,z= 4970000)
BBOX <- 50000
GPS_UTC_OFFSET <- as.integer(as.POSIXct('1980-01-06',tz="UTC"))
GPS_UTC_LEAPSECONDS <- -18

nanosToDate <- function(time){
  return(as.Date(as.POSIXct(time/1e9+GPS_UTC_OFFSET+GPS_UTC_LEAPSECONDS,origin = "1970-01-01", tz = "UTC")))
}

UTCtoNanos <- function(UTCTime) {
  return((as.numeric(UTCTime) -GPS_UTC_OFFSET -GPS_UTC_LEAPSECONDS)*1e9)
}

makeSpatial <- function(gnssData){
  zeroMissingAltitudes <- function(epochs){
    missing <- which(is.nan(epochs[["Altitude"]]))
    if (length(missing)!=0) {
      warning(paste0(length(missing)," altitudes were missing, and set to 0",collapse=""))
      epochs[missing,"Altitude"]=0
    }
    return(epochs)
  }
  
  measurementsSF <- locateSatellites(gnssData)
  epochsSF <- gnssData[["epochs"]] %>% zeroMissingAltitudes() %>% st_as_sf(coords=c("Longitude","Latitude","Altitude"),dim="XYZ",crs=CRS_LONGLAT_WGS84) %>% st_transform(CRS_CART_WGS84)
  index <- match(measurementsSF[["EpochID"]],epochsSF[["EpochID"]])
  delta <- st_coordinates(measurementsSF)-st_coordinates(epochsSF)[index,]
  measurementsSF[["SatDistance"]]  <- apply(delta^2,1,sum)^0.5
  measurementsBNG <- measurementsAsBNG(measurementsSF)
  epochsBNG <- epochsSF %>% st_transform(CRS_BNG)
  return(list(measurements=measurementsBNG,epochs=epochsBNG))
}

locateSatellites <- function(gnssData){
  #corrects for satellite clock errors and adds satellite locations to a set of measurements.The epoch datatable is also required as it provides epoch times.
  
  measurements <- gnssData[["measurements"]]
  epochs <- gnssData[["epochs"]]
  locationFunctions <- as.Date(epochs[["UTCTime"]]) %>% unique() %>% set_names() %>% map(getSVIDLocation) #This creates the getSVIDLocation(date) functions in this environment
  
  isEstimateMeasurementTime <- is.na(measurements[["TransmitterGPSTimeNanos"]])
  imputedReceiverGPSTimeNanos <- ifelse(is.na(measurements[["ReceiverGPSTimeNanos"]]),
                                          epochs[["UTCTime"]][match(measurements[["EpochID"]],epochs[["EpochID"]])] %>% UTCtoNanos(),
                                      measurements[["ReceiverGPSTimeNanos"]])
  
  satelliteLocations <- pmap(list(imputedReceiverGPSTimeNanos %>% nanosToDate() %>% as.character(), isEstimateMeasurementTime, measurements[["TransmitterGPSTimeNanos"]],imputedReceiverGPSTimeNanos ,measurements[["Svid"]]),
                             ~locationFunctions[[..1]](..2,..3,..4,..5)) %>% transpose() %>% simplify_all %>% as_tibble()
  satelliteLocations <- satelliteLocations %>% mutate(rotation=(Rx^2+Ry^2+Rz^2)^0.5,velocity=(Vx^2+Vy^2+Vz^2)^0.5)

  measurementsSF <- bind_cols(measurements,satelliteLocations) %>% filter(!is.na(x)) %>%
    st_as_sf(coords=c("x","y","z"),dim="XYZ",crs=CRS_CART_WGS84)
  cat(nrow(measurementsSF),"satellites located out of a potential",nrow(measurements),"\n",
      sum(!is.na(measurementsSF[["ReceivedSvTimeNanos"]]))," received satellites located out of a potential",sum(!is.na(measurements[["ReceivedSvTimeNanos"]])))
  return(measurementsSF)
}

measurementsAsBNG <- function(measurements) {
  #This bounds the measurements to a 100km box around London, to ensure validity of the transformationm, before transforming to BNG. THe 27700 transformation is accurate to within a few metres. This is acceptable as the error is scaled down proportionally when consider the intesection of a device and building)
  clip_to_box <- function(point){
    delta <- point-ORIGIN
    if (max(abs(delta))<1e-8) {
      return(ORIGIN)
    } else {
      scale <- BBOX/max(abs(delta))
      return(ORIGIN + scale*delta )
    }
  }
  
  coords <- st_coordinates(measurements)
  new_coords <- map(array_branch(coords,1),clip_to_box) %>% transpose %>% simplify_all() %>% as_tibble
  st_geometry(measurements) <- st_as_sf(new_coords,coords=c("x","y","z"),dim="XYZ",crs=CRS_CART_WGS84) %>% st_geometry()
  measurements <- st_transform(measurements,CRS_BNG)
  return(measurements)
}

combineAsLines <- function(gnssData){
  #merges the measurement and epochs data, turning device and satellite points into line geometry, and adding device-centred azimuth and elevation info)
  coords_to_line <- function(x,y){
    return(paste0("LINESTRING Z (",x[1],
                  " ",x[2],
                  " ",x[3],
                  ", ",y[1],
                  " ",y[2],
                  " ",y[3],
                  ")")
    )
  }
  
  measurements <- gnssData[["measurements"]]
  epochs <- gnssData[["epochs"]]
  crs <- st_crs(measurements)
  if(st_crs(epochs)!=crs){
    message("CRS conversion was required")
    epochs <- st_transform(epochs,crs)
  }
  
  epochsToJoin <- epochs[match(measurements[["EpochID"]],epochs[["EpochID"]]),]
  relativeDirection <- getAzimuthELevation(epochsToJoin[["geometry"]],measurements[["geometry"]])
  satelliteCoords <- st_coordinates(measurements)
  deviceCoords <- st_coordinates(epochsToJoin)
  lines <- map2(array_branch(deviceCoords,1),array_branch(satelliteCoords,1),~coords_to_line(.x,.y)) %>% st_as_sfc(crs=crs)
  df <- bind_cols(st_set_geometry(measurements,lines),st_drop_geometry(epochsToJoin[,"UTCTime"]),relativeDirection)
  return(df)
}

addAtmosphericDelay <- function(lines){
  isModelled <- !is.na(lines[["TransmitterGPSTimeNanos"]])
  LLA <- lines %>% subset(isModelled) %>% st_transform(CRS_LONGLAT_WGS84) %>% st_coordinates() %>% .[seq.int(1,nrow(.),by=2),-4]
  colnames(LLA) <- c("long","lat","altitude")
  times <- lines[["TransmitterGPSTimeNanos"]] %>% subset(isModelled)
  dates <-  times%>% nanosToDate()%>% as.character()
  elevations <- lines[["elevation"]]%>% subset(isModelled) 
  azimuths <- lines[["azimuth"]] %>% subset(isModelled) 
  ionosphereFunctions <- dates %>% unique() %>% set_names() %>% map(getIonosphericDelay) #creates the functions in this environment
  troposphereDelay=getTroposphericDelay()
  
  lines[["ionosphereDelay"]] <- NA_real_
  lines[isModelled,"ionosphereDelay"] <- pmap_dbl(list(dates,elevations,azimuths,LLA[,"lat"],LLA[,"long"],times),
                                                  ~ionosphereFunctions[[..1]](elevation= ..2, azimuth= ..3,latitude= ..4,longitude= ..5,GPSTimeNanos= ..6))
  lines[["troposphereDelay"]] <- NA_real_
  lines[isModelled,"troposphereDelay"] <- pmap_dbl(list(elevations,LLA[,"altitude"],LLA[,"lat"],dates),
                                                   ~troposphereDelay(elevation = ..1,height = ..2,latitude = ..3,date = ..4))

  return(lines)  
}

getAzimuthELevation <- function(device,satellite){
  #reports wrt to device 
  deviceENU <- st_transform(device,CRS_BNG)%>% st_coordinates()
  satENU <-  st_transform(satellite,CRS_BNG)%>% st_coordinates() 
  delta <- (satENU-deviceENU) %>% as_tibble()
  azel <-delta %>% transmute(elevation=asin(Z/(X^2+Y^2+Z^2)^0.5)*180/pi,azimuth=atan2(X,Y)*180/pi)
  return(azel)
}

intersectBuilding <- function(polygonSF,linesSF){
  #this takes the dataset of satellite observations and intersects them with the building (the usual sf function does not account for height properly or ensure the closest point to the origin is taken).
  #It also calculates 4 useful parameters for analysis of data spatial dependencies x,y,z,w:
  
  intersectPoints <- function(lines,polygon){
    returnPartLineXYZ <- function (lineXYZ,partLineXY){
      z <- lineXYZ[1,3] + (lineXYZ[2,3]-lineXYZ[1,3])/(lineXYZ[2,1]-lineXYZ[1,1])*(partLineXY[2,1]-partLineXY[1,1])
      return(c(partLineXY[2,1],partLineXY[2,2],z))
    }
    
    intersectingPoints <- list(st_point(c(NA_real_,NA_real_,NA_real_),dim="XYZ"))[rep(1,NROW(lines))] %>% st_as_sfc() 
    isIntersecting <- st_intersects(lines,polygon,sparse=FALSE) %>% drop()
    if (all(!isIntersecting)) {return(NULL)}
    originalCoordsXYZ <- lines[isIntersecting] %>% st_coordinates() %>% split.data.frame(.[,4])
    intersectingMultiLinestrings <- lines[isIntersecting] %>% st_difference(polygon) %>% st_cast("MULTILINESTRING")
    intersectingCoordsXY <- intersectingMultiLinestrings %>% st_coordinates() %>% .[.[,3]==1,]%>% split.data.frame(.[,4])
    intersectingPoints[isIntersecting]  <- map2(originalCoordsXYZ,intersectingCoordsXY,~returnPartLineXYZ(.x,.y)) %>% map(st_point,dim="XYZ")
    a <- intersectingPoints%>% st_set_crs(st_crs(lines))
    return (a)
  }

  getDistanceAlong <- function(polygon,points){
    isIntersecting <- map_lgl(points,~!st_is_empty(.x))
    distances <- as.double(rep(NA_real_,times=NROW(points)))
    bufferedPoints <- points[isIntersecting] %>% st_buffer(1e-2)
    distances[isIntersecting] <- st_difference(polygon,bufferedPoints) %>% map(~.x[[1]] %>% st_linestring()) %>% map_dbl(st_length)
    return(distances)
  }
  polygon <- st_geometry(polygonSF) 
  lines <- st_geometry(linesSF)
  linesSF[["intersectingPoints"]] <- intersectPoints(lines,polygon)
  linesSF[["X"]] <- map_dbl(linesSF[["intersectingPoints"]],~.x["X"])
  linesSF[["Y"]] <- map_dbl(linesSF[["intersectingPoints"]],~.x["Y"])
  linesSF[["Z"]] <- map_dbl(linesSF[["intersectingPoints"]],~.x["Z"])
  # linesSF[["W"]] <- getDistanceAlong(polygon %>% st_cast("LINESTRING"),linesSF[["intersectingPoints"]])
  linesSF[["X0"]] <- map_dbl(linesSF[["geometry"]],~.x[[1]])
  linesSF[["Y0"]] <- map_dbl(linesSF[["geometry"]],~.x[[3]])
  linesSF[["Z0"]] <- map_dbl(linesSF[["geometry"]],~.x[[5]])
  linesSF[["X1"]] <- map_dbl(linesSF[["geometry"]],~.x[[2]])
  linesSF[["Y1"]] <- map_dbl(linesSF[["geometry"]],~.x[[4]])
  linesSF[["Z1"]] <- map_dbl(linesSF[["geometry"]],~.x[[6]])
  return(linesSF)
}

