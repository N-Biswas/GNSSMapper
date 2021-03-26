require(tidyverse)
source("code/processing/logTools.R")
source("code/processing/sfTools.R")


# 1. logTools:reads logger data into a datatable of epochs (device location and time of measurement) and measurements (gnss raw measurements by sattelite and epoch). 
# 2. sfTools: converts epochs into a sf points dataframe
#             converts measurements into a sf points dataframe. Uses epoch times from epoch database and sp3 functions which obtains ephemeris information from ESA
#             amalgamates epochs and measurements sf into lines from device to satellite and creates a bounding box
# 3. have map on UCL and intersect to find interesting measurements

#TO-DO
# create a different set of coords for line data (X,Y,Z,W)
# calculate precise GNSS pseudorange. add GNSS state to avoid missing data. There are sat dependent errors - thinking this is down to gps time versus utc time and clock error at the satellite level (how is clock error identified in broadcasts?)
#I have the satellite precise location and the time of epoch :
# 1) I use transmitter stamp to work backwards to transmit time  when state is high enough. I also need a sat position for low/no states, but this is estimated based on sat altitude.
# 2) timestamps corrected using epoch data
# 3) # correct for relativity term in clock bias
# correcot for atmospheric effects

#correct for sagnac effect

#currently dropping glonass satellites which return a frequency channel number. Need to convert them into a OCN.
#why can some received satellites not be found in the ephemeris? I assume non received is due to their being out of service
# build a system to turn classifiers into a point

processFile <- function(file,overwrite=FALSE,outputFolder,comment=NA,calculatePR=FALSE,observerLocationMethod=1,LocationList=NULL){
  #takes a logger file and produces a dataframe of recorded signals. arguments are:
  #overwrite existing file?
  #which folder to put the prcoessed file (same name as logger file)
  #comment to add to related metadata file
  #calculate psuedorange or just signal strength
  #method to give observer lcoation 1 - use LocationList, 2 -use GPS fixes but use same epochs as the location list, 3 -use GPS fixes are keep all possible epochs
  
  filepath <- paste0(getwd(),"/data/",file,collapse="")
  filename <- substr(file,1,nchar(file)-4)
  folder <- paste0(getwd(),"/processedData/",outputFolder,"/",collapse="")
  writepath <- paste0(folder,filename,".csv.gz",collapse="")
  metapath <-  paste0(folder,filename,"_meta.csv",collapse="")
  
  if (!dir.exists(folder)) {dir.create(folder)}
  
  if(file.exists(metapath) & !overwrite) {return(NULL)}
  
  loggerOutput <- switch(observerLocationMethod,
                         readLog(filepath,calculatePR=calculatePR,onlyMarkedLocations=TRUE,LocationList=LocationList,useFixes=FALSE),
                         readLog(filepath,calculatePR=calculatePR,onlyMarkedLocations=TRUE,LocationList=LocationList,useFixes=TRUE),
                         readLog(filepath,calculatePR=calculatePR,onlyMarkedLocations=FALSE,LocationList=LocationList,useFixes=TRUE))
  
  gnssData <- loggerOutput[1:2]
  metaData <- bind_cols(tibble(comment=comment,calculatePR=calculatePR,observerLocationMethod=observerLocationMethod),loggerOutput[[3]])
  
  if (NROW(gnssData$measurements)==0) {
    write_csv(metaData,metapath)
    return(NULL)
  }
  
  gnssSF <- makeSpatial(gnssData)
  gnssLines <- combineAsLines(gnssSF)
  groundFilter <- gnssLines$elevation<0
  vertFilter <- gnssLines$elevation>85
  aboveGround <- gnssLines[!groundFilter & !vertFilter,]
  if(calculatePR){aboveGround <- addAtmosphericDelay(aboveGround)}
  output <- intersectBuilding(UCL,aboveGround)
  output <- output %>% st_drop_geometry() %>% select(-intersectingPoints)
  metaData <- metaData %>% mutate(receivedSignals=sum(!is.na(output[["ReceivedSvTimeNanos"]])),missingSignals=sum(is.na(output[["ReceivedSvTimeNanos"]])))
  
  write_csv(output,writepath)
  write_csv(metaData,metapath)
}


MAP_PATH <- "/map/"
MAP_FILE <- "mastermap-topo_3473984_0.gml.gz"
LAYER <- "TopographicArea"
filename <- paste0(getwd(),MAP_PATH,MAP_FILE,collapse = "")
mastermap <- st_read(filename,LAYER)
# bbox <- st_bbox(c(xmin=529450,xmax=529700,ymin=182150,ymax=182420))
# boundedMastermap <- mastermap%>% st_crop(bbox)
UCL <- mastermap %>% filter(fid=="osgb5000005153747520")
rm(mastermap)



files <- list.files(path=paste0(getwd(),"/data",collapse=""),
                    pattern="*.txt",
                    full.name=FALSE,
                    recursive=TRUE)
ZERO_HEIGHT <- 80
LocationList <- read_csv(paste0(getwd(),"/data/","location_labels_20_02_13.csv",collapse="")) %>% transmute(LocationID=LocationID,Latitude=lat,Longitude=long,Altitude=ZERO_HEIGHT)

for (file in files) {
  processFile(file,overwrite=TRUE,outputFolder="gnss_exact_location",LocationList = LocationList)
}

for (file in files) {
  processFile(file,overwrite=TRUE,outputFolder="gnss_estimated_location_small",LocationList = LocationList,observerLocationMethod = 2)
}

for (file in files) {
  processFile(file,overwrite=TRUE,outputFolder="gnss_estimated_location_all",LocationList = LocationList,observerLocationMethod = 3)
}


joinProcessedFiles <- function(folder){
  processedFiles <- list.files(path=paste0(getwd(),"/processedData/",folder,collapse=""),
                               pattern="*.csv.gz",
                               full.name=TRUE,
                               recursive=TRUE)
  metaFiles <- list.files(path=paste0(getwd(),"/processedData/",folder,collapse=""),
                               pattern="*_meta.csv",
                               full.name=TRUE,
                               recursive=TRUE)
  
  gnss <- map(processedFiles,read_csv) %>% data.table::rbindlist()
  metaData <- map(metaFiles,read_csv) %>% data.table::rbindlist(fill=TRUE)
  short <- gnss %>% select(EpochID,Svid,Cn0DbHz,UTCTime,elevation,X,Y,Z,X0,Y0,Z0,X1,Y1,Z1)
  
  write_csv(short,paste0(getwd(),"/processedData/",folder,".csv.gz",collapse=""))
  write_csv(metaData,paste0(getwd(),"/processedData/",folder,"_meta.csv",collapse=""))
}

joinProcessedFiles("gnss_exact_location")
joinProcessedFiles("gnss_estimated_location_small")
joinProcessedFiles("gnss_estimated_location_all")




