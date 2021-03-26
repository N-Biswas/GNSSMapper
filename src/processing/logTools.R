require(tidyverse)

#tools for extracting information from the logger output files
# 1. A epoch datatable - a series of epochs with location of device and time. 
# 2. a measurements datatable for each epoch and all gnss sats:e.g. pseudorange, doppler and signalstrength



readLog <- function(file,calculatePR,...){
    cat("Reading",file,"\n")
    data <- readGnss(file) %>% preprocessGnss(calculatePR)
    epochs <- getEpochs(data,...) 
    NumberAllEpochs <- NROW(epochs)
    numberAllSignals <- NROW(data[["Raw"]])
    epochs <- epochs[!missing(epochs) & !is.na(epochs$UTCTime),] #ensures valid time and position measurement
    measurements <- data[["Raw"]] %>% filter(EpochID %in% epochs[["EpochID"]]) %>% getMeasurements()
    metaData <- tibble(AllEpochs=NumberAllEpochs,AllSignals=numberAllSignals,RecordedEpochs=NROW(epochs))
    gnssData <- list(measurements=measurements,epochs=epochs,metaData=metaData)
  return(gnssData)
}

stateCheck <- function(state,number) {
  asReverseBinary <- function(x){
    if (x<=1){
      return (x %% 2)
    } else {
      return(c(x %% 2,asReverseBinary  (x %/% 2)))
    }
  }
  
  bin <- asReverseBinary(state)
  test <- length(asReverseBinary(number))
  return(length(bin)>=test && bin[test]==1)  
  
}
readGnss <- function(filename){
  #Reads a log file and returns 2 dataframes: raw signals, and associated epoch)
  DATA_TYPES <- c(Raw="Raw",Fix="Fix")
  HEADERS<- list(Raw=c("ElapsedRealtimeNanos","TimeNanos","FullBiasNanos","BiasNanos","BiasUncertaintyNanos","Svid","TimeOffsetNanos","State","ReceivedSvTimeNanos","ReceivedSvTimeUncertaintyNanos","Cn0DbHz","PseudorangeRateMetersPerSecond","Constellation","UserLatitude","UserLongitude","UserAltitude"),
                 Fix=c("Latitude","Longitude","Altitude","UTCTimeInMs","ElapsedRealtimeNanos"))
  COLUMNS <- list(Raw=c(2,3,6,7,8,12,13,14,15,16,17,18,29,32,33,34),Fix=c(3,4,5,8,9))
  
  extractData <- function(data,type){
    a <- data[startsWith(data,type)]
    if(length(a)==0) {
      return(NULL)
    } else {
      b <- read_csv(file=a,col_names = FALSE)
      c <- b[COLUMNS[[type]]]
      colnames(c) <- HEADERS[[type]]
      return(c)
    }
  } 
  
  log <- readLines(filename)
  data <- map(DATA_TYPES,~extractData(log,.))
  return(data)
}

calculateTransmitterTimeNanos <- function(ReceivedSvTimeNanos,ConstellationLetter,State,ReceiverGPSTimeNanos){
  calculateGPSTravelTime <- function(ReceivedSvTimeNanos,State,ReceiverGPSTimeNanos) {
    if(!stateCheck(State,8)) {return(NA_real_)} #requiring full TOW with no measurement ambiguity. 
    NUMBER_NANOS_WEEK <-604800*10^9 
    ReceiverTimeNanos <- ReceiverGPSTimeNanos - NUMBER_NANOS_WEEK * (ReceiverGPSTimeNanos %/% NUMBER_NANOS_WEEK)
    
    return (ReceiverGPSTimeNanos-(ReceiverTimeNanos-ReceivedSvTimeNanos))
  }
 
  calculateGalileoTravelTime <- function(ReceivedSvTimeNanos,State,ReceiverGPSTimeNanos) {
    if(!stateCheck(State,8)| !stateCheck(State,2048) ) {return(NA_real_)} #requiring second code lock as well as TOW decoded.
    NUMBER_NANOS_100MILLIS <-10^8 
    ReceiverTimeNanos <- ReceiverGPSTimeNanos - NUMBER_NANOS_100MILLIS * (ReceiverGPSTimeNanos %/% NUMBER_NANOS_100MILLIS)
    ReceivedSvPilotTimeNanos <- ReceivedSvTimeNanos - NUMBER_NANOS_100MILLIS * (ReceivedSvTimeNanos %/% NUMBER_NANOS_100MILLIS) #Documentation suggests that the measurement can collapse into measuring pilot stage with 100ms ambiguity. This step ensures conformity of measurment
    return (ReceiverGPSTimeNanos-(ReceiverTimeNanos-ReceivedSvPilotTimeNanos))
  } 
  
  calculateBeidouTravelTime <- function(ReceivedSvTimeNanos,State,ReceiverGPSTimeNanos) {
    if(!stateCheck(State,8)) {return(NA_real_)} #requiring full TOW with no measurement ambiguity. 
    NUMBER_NANOS_WEEK <-604800*10^9 
    BEIDOU_OFFSET <- 14*10^9
    ReceiverTimeNanos <- ReceiverGPSTimeNanos - NUMBER_NANOS_WEEK * (ReceiverGPSTimeNanos %/% NUMBER_NANOS_WEEK) - BEIDOU_OFFSET
    return (ReceiverGPSTimeNanos-(ReceiverTimeNanos-ReceivedSvTimeNanos))
  }
  
  calculateGlonassTravelTime <- function(ReceivedSvTimeNanos,State,ReceiverGPSTimeNanos) {
    if(!stateCheck(State,128)) {return(NA_real_)} #requiring full TOD with no measurement ambiguity. 
    NUMBER_NANOS_DAY <-86400*10^9 
    GLONASS_3HOUR_OFFSET <- 10800*10^9
    GLONASS_LEAPSECOND_OFFSET <- -18*10^9
    ReceiverTimeNanos <- ReceiverGPSTimeNanos - NUMBER_NANOS_DAY * (ReceiverGPSTimeNanos %/% NUMBER_NANOS_DAY) + GLONASS_3HOUR_OFFSET + GLONASS_LEAPSECOND_OFFSET
    return (ReceiverGPSTimeNanos-(ReceiverTimeNanos-ReceivedSvTimeNanos))
  }
  
  TransmitterGPSTimeNanos <- pmap_dbl(list(ConstellationLetter,ReceivedSvTimeNanos,State,ReceiverGPSTimeNanos),~switch(..1,
                                                                                                          "G"=calculateGPSTravelTime(..2,..3,..4),
                                                                                                          "R"=calculateGlonassTravelTime(..2,..3,..4),
                                                                                                          "C"=calculateBeidouTravelTime(..2,..3,..4),
                                                                                                          "E"=calculateGalileoTravelTime(..2,..3,..4),
                                                                                                          NA_real_))
  return(TransmitterGPSTimeNanos)
}

preprocessGnss <- function(rawData,calculatePR){
  CONSTELLATION <- c("G"=1,"R"=3,"C"=5,"E"=6) #GPS,GLONASS,BEIDOU,GALILEO
  GPS_UTC_OFFSET <- as.integer(as.POSIXct('1980-01-06',tz="UTC"))
  GPS_UTC_LEAPSECONDS <- -18
  epochTimes <- unique(rawData[["Raw"]][["ElapsedRealtimeNanos"]]/1e9)
  
  Raw <- 
    rawData[["Raw"]] %>% 
    mutate(Svid = paste0(names(CONSTELLATION)[match(Constellation,CONSTELLATION)],formatC(Svid,width=2,flag="0")),
           UTCTime = as.POSIXct(as.integer((TimeNanos-(FullBiasNanos+BiasNanos))/1e9) +GPS_UTC_OFFSET+GPS_UTC_LEAPSECONDS,origin = "1970-01-01", tz = "UTC"),
           ReceiverGPSTimeNanos = ifelse(FullBiasNanos=="",NA_real_,TimeNanos-(FullBiasNanos+BiasNanos)+TimeOffsetNanos), #This checks for a valid time fix for the epoch
           TransmitterGPSTimeNanos = ifelse(calculatePR,calculateTransmitterTimeNanos(ReceivedSvTimeNanos,names(CONSTELLATION)[match(Constellation,CONSTELLATION)],State,ReceiverGPSTimeNanos),NA_real_), #this converts the timestamp of the sv transmission into a GPS time format, if ambiguities are resolvable
           SystemTime=ElapsedRealtimeNanos/1e9,
           EpochID=match(SystemTime,epochTimes),
           UserLongitude=ifelse(is.nan(UserLongitude),NA,UserLongitude),
           UserAltitude=ifelse(is.nan(UserAltitude),NA,UserAltitude)) %>% 
    select(-TimeNanos,-TimeOffsetNanos,-ElapsedRealtimeNanos,-FullBiasNanos,-BiasNanos,-Constellation) 
  
  if(is_null(rawData[["Fix"]])){
    Fix <- rawData[["Fix"]]
  } else {
    Fix <- 
      rawData[["Fix"]] %>% 
      mutate(UTCTime = as.POSIXct(UTCTimeInMs/1e3,origin = "1970-01-01", tz = "UTC"),
             SystemTime=ElapsedRealtimeNanos/1e9) %>%
      select(SystemTime,UTCTime,Latitude,Longitude,Altitude) 
  }  
  return(list(Raw=Raw,Fix=Fix))
}

# test <- preprocessGnss( readGnss(files[1]),0)
# raw <- test[["Raw"]]
# summary(raw)
# table(raw %>% filter(!is.na(PseudoRangeMetres)) %>% select(Constellation))
# plot()
# gps <- raw %>% filter(Constellation==1)
# summary(raw %>% filter(Constellation!=3))
# raw[["Bin"]]<- intToBin(raw$State) 
# plot(raw$Constellation,raw$PseudoRangeNanos)
# gps[["ReceiverTimeOFWeeksNanos"]]=gps$GPSTimeNanos- NUMBER_NANOS_WEEK *(gps$GPSTimeNanos %/% NUMBER_NANOS_WEEK )
# gps[[]]
# summary <- gps %>% select(ReceiverTimeOFWeeksNanos,ReceivedSvTimeNanos)
# summary(gps)

getMeasurements <- function(rawData){
  epochIDs <- unique(rawData[["EpochID"]])
  svids <- c(c("E01","E02","E03","E04","E05","E06","E07","E08","E09","E10","E11","E12","E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25","E26","E27","E28","E29","E30","E31","E32","E33","E34","E35","E36"),
            c("G01","G02","G03","G04","G05","G06","G07","G08","G09","G10","G11","G12","G13","G14","G15","G16","G17","G18","G19","G20","G21","G22","G23","G24","G25","G26","G27","G28","G29","G30","G31","G32"),                   
            c("R01","R02","R03","R04","R05","R06","R07","R08","R09","R10","R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","R21","R22","R23","R24"),                      
            c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31","C32","C33","C34","C35","C36","C37"))
  measurements <- tibble(EpochID = rep(epochIDs,each=length(svids)),
                   Svid = rep(svids,times=length(epochIDs)))
  
  measurements <- left_join(measurements,rawData,by=c("EpochID","Svid")) %>% select(-UserLatitude,-UserLongitude,-UserAltitude,-SystemTime,-UTCTime)
  drop <- nrow(rawData)-sum(!is.na(measurements[["ReceivedSvTimeNanos"]]))
  
  cat(sum(!is.na(measurements[["ReceivedSvTimeNanos"]])),"signals received out of a potential",nrow(measurements),"\n",
      drop," additional signals were received but could not be matched to our svid database (covers GPS,Galileo,Glonass and Beidou)\n"
  )
  
  return(measurements)
}

getEpochs <- function(gnssData,onlyMarkedLocations=TRUE,LocationList=NULL,useFixes=FALSE){
  #This function produces a data table of time and location for each epoch. If the location is provided by the user it is used. 
  #userLocationList indicates whether the LLA has been fully specified or whether a location number has been provided (in Latitude column) to be looked up iun a list.
  #if usesFixes is true, any epochs without a user-added location will be estimated based on the GNSS reciever measurements.
  #Any epochs without a location are dropped
  epochs <- distinct(gnssData[["Raw"]],EpochID, .keep_all=TRUE) %>% select(EpochID,SystemTime,UTCTime,Latitude=UserLatitude,Longitude=UserLongitude,Altitude=UserAltitude)
  epochs <- writeUserLLA(epochs,LocationList)
  if (useFixes){
    epochs[!missing(epochs),] <- getDeviceFixes(gnssData[["Fix"]],epochs[!missing(epochs),])
    if(!onlyMarkedLocations){epochs[missing(epochs),] <- getDeviceFixes(gnssData[["Fix"]],epochs[missing(epochs),])}
    
  }
  
return(epochs)
}

missing <- function(epochs){
  return((is.na(epochs[["Latitude"]]) | is.na(epochs[["Longitude"]]) | is.na(epochs[["Altitude"]])))
}

getDeviceFixes <- function(fix,epochs){
  fix <- fix[order(fix[["UTCTime"]]),] #required in order to use findInterval
  fixByEpochIndex  <- findInterval(epochs[["UTCTime"]],fix[["UTCTime"]],all.inside = TRUE)
  fixByEpoch <- fix[fixByEpochIndex,] 
  epochs <- cbind(epochs %>% select(EpochID),fixByEpoch)
  return(epochs)
}

writeUserLLA <- function(epochs,LocationList){
  epochs <- epochs %>% mutate(LocationID=Latitude) %>% select(-Latitude,-Longitude,-Altitude)
  epochs <- left_join(epochs,LocationList,by="LocationID") %>% select(-LocationID)
  return(epochs)
}

# addToDatabase <- function(gnssData,gnssDatabase=NULL){
#   if(is.null(gnssDatabase)) {
#     epochOffset <- 0
#   } else{
#     epochOffset <- max(gnssDatabase[["epochs"]][["EpochID"]])
#   }
#   gnssData[["measurements"]][["EpochID"]] <- gnssData[["measurements"]][["EpochID"]]+epochOffset
#   gnssData[["epochs"]][["EpochID"]] <- gnssData[["epochs"]][["EpochID"]]+epochOffset 
#   gnssDatabase[["measurements"]] <- bind_rows(gnssDatabase[["measurements"]],gnssData[["measurements"]])
#   gnssDatabase[["epochs"]] <- bind_rows(gnssDatabase[["epochs"]],gnssData[["epochs"]])
#   return(gnssDatabase)
# }


