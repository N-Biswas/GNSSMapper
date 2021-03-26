require(tidyverse)
require(httr)

RINEX_LOCAL <- "ephemeris"
USERNAME <- ""
PASSWORD <- ""

dateToYearDay <- function(input){
  date<- as.Date(input)
  year <- format(date,"%Y")
  YearDay <- c(year = format(date,"%Y"),
               day = format(date,"%j"))
  return(YearDay)
}

getRinexFile <- function(date){
  RINEX_DATASITE <- "https://cddis.nasa.gov/archive/gnss/data/daily/"
  rinexDate <- dateToYearDay(date)
  filename <- paste0("BRDC00IGS_R_",rinexDate["year"],rinexDate["day"],"0000_01D_MN.rnx.gz",collapse="")
  url <- paste0(c(RINEX_DATASITE,rinexDate["year"],"/brdc/",filename),collapse="")
  local <- paste0(c(getwd(),"/",RINEX_LOCAL,"/",filename),collapse="")
  
  cat(url)
  if (!dir.exists(RINEX_LOCAL)) {dir.create(RINEX_LOCAL)}
  if (!file.exists(local)) {
    GET(url,write_disk(local),authenticate(USERNAME,PASSWORD))
  }
  
  conn <- gzfile(local)
  rinex <- readLines(conn)
  close(conn)
  return(rinex)
}


getNavigation <- function(rinex){
  #from a rinex file with format here: ftp://igs.org/pub/data/format/rinex303.pdf
  #return a tibble with navigation parameters by time and svid
  getMessage <- function(x){
    stringToUTC <- function(string,constellation){
      return(string)
    }
    
    readGPS <- function(x){
      return(
        list(Svid = str_sub(x,1,3),
             UTCtime = stringToUTC(str_sub(x,5,23),"GPS"),
             Kepler = as.Kepler(timeOfEphemeris = as.numeric(str_sub(x,245,263)),
                                sqrtMajorAxis = as.numeric(str_sub(x,222,240)),
                                eccentricity = as.numeric(str_sub(x,184,202)),
                                meanAnomaly0 = as.numeric(str_sub(x,142,160)),
                                argumentOfPerigree  = as.numeric(str_sub(x,363,381)),
                                inclination0 = as.numeric(str_sub(x,325,343)),
                                omega0 = as.numeric(str_sub(x,283,301)),
                                meanMotionDifference  = as.numeric(str_sub(x,123,141)),
                                inclinationRate = as.numeric(str_sub(x,405,423)),
                                omegaRate = as.numeric(str_sub(x,382,400)),
                                cuc = as.numeric(str_sub(x,165,183)),
                                cus = as.numeric(str_sub(x,203,221)),
                                crc = as.numeric(str_sub(x,344,362)),
                                crs = as.numeric(str_sub(x,104,122)),
                                cic = as.numeric(str_sub(x,264,282)),
                                cis = as.numeric(str_sub(x,302,320)),
                                clockOffset = as.numeric(str_sub(x,24,42)),
                                clockDrift = as.numeric(str_sub(x,43,61)),
                                clockDriftRate = as.numeric(str_sub(x,62,80))),
             Cartesian = NA,
             # IODE = as.numeric(str_sub(x,85,103)),
             # Codes = as.numeric(str_sub(x,424,442)),
             # GPSWeek = as.numeric(str_sub(x,443,461)),
             # L2P_flag = as.numeric(str_sub(x,462,480)),
             # Accuracy = as.numeric(str_sub(x,485,503)),
             # Health = as.numeric(str_sub(x,504,522)),
             # TGD = as.numeric(str_sub(x,523,541)),
             # IODC = as.numeric(str_sub(x,542,560)),
             # TransmissionTime = as.numeric(str_sub(x,565,583)),
             # FitInterval = as.numeric(str_sub(x,584,602))
        )
      )
    }
    
    readGalileo <- function(x){
      return(
        list(Svid = str_sub(x,1,3),
             UTCtime = stringToUTC(str_sub(x,5,23),"Gal"),
             Kepler = as.Kepler(timeOfEphemeris = as.numeric(str_sub(x,245,263)),
                                sqrtMajorAxis = as.numeric(str_sub(x,222,240)),
                                eccentricity = as.numeric(str_sub(x,184,202)),
                                meanAnomaly0 = as.numeric(str_sub(x,142,160)),
                                argumentOfPerigree  = as.numeric(str_sub(x,363,381)),
                                inclination0 = as.numeric(str_sub(x,325,343)),
                                omega0 = as.numeric(str_sub(x,283,301)),
                                meanMotionDifference  = as.numeric(str_sub(x,123,141)),
                                inclinationRate = as.numeric(str_sub(x,405,423)),
                                omegaRate = as.numeric(str_sub(x,382,400)),
                                cuc = as.numeric(str_sub(x,165,183)),
                                cus = as.numeric(str_sub(x,203,221)),
                                crc = as.numeric(str_sub(x,344,362)),
                                crs = as.numeric(str_sub(x,104,122)),
                                cic = as.numeric(str_sub(x,264,282)),
                                cis = as.numeric(str_sub(x,302,320)),
                                clockOffset = as.numeric(str_sub(x,24,42)),
                                clockDrift = as.numeric(str_sub(x,43,61)),
                                clockDriftRate = as.numeric(str_sub(x,62,80))),
             Cartesian = NA,
             # IODnav = as.numeric(str_sub(x,85,103)),
             # DataSources = as.numeric(str_sub(x,424,442)), #float to int
             # GalWeek = as.numeric(str_sub(x,443,461)),
             # # L2P_flag = as.numeric(str_sub(x,462,480)),
             # Accuracy = as.numeric(str_sub(x,485,503)),
             # Health = as.numeric(str_sub(x,504,522)), #float to int
             # BGD_E5a = as.numeric(str_sub(x,523,541)),
             # BGD_E5b = as.numeric(str_sub(x,542,560)),
             # TransmissionTime = as.numeric(str_sub(x,565,583))
        )
      )
    }
    
    readGlonass <- function(x){
      return(
        list(Svid = str_sub(x,1,3),
             UTCtime = stringToUTC(str_sub(x,5,23),"Glo"),
             Kepler = NA,
             Cartesian = as.Cartesian(timeOfEphemeris,
                                      x = as.numeric(str_sub(x,85,103)),
                                      y = as.numeric(str_sub(x,165,183)),
                                      z = as.numeric(str_sub(x,245,263)),
                                      v_x = as.numeric(str_sub(x,104,122)),
                                      v_y = as.numeric(str_sub(x,184,202)),
                                      v_z = as.numeric(str_sub(x,264,282)),
                                      a_x = as.numeric(str_sub(x,123,141)),
                                      a_y = as.numeric(str_sub(x,203,221)),
                                      a_z = as.numeric(str_sub(x,283,301)),
                                      clockOffset =  as.numeric(str_sub(x,24,42)),
                                      clockFrequencyOffset = as.numeric(str_sub(x,43,61))),
             # MessageFrameTime = as.numeric(str_sub(x,62,80)),
             # Health = as.numeric(str_sub(x,142,160)),
             # Frequency_number = as.numeric(str_sub(x,222,240)),
             # Age_Info = as.numeric(str_sub(x,302,320))
        )     
      )
    }
    
    readBeidou <- function(x){
      return(
        list(Svid = str_sub(x,1,3),
             UTCtime = stringToUTC(str_sub(x,5,23),"Bei"),
             Kepler = as.Kepler(timeOfEphemeris = as.numeric(str_sub(x,245,263)),
                                sqrtMajorAxis = as.numeric(str_sub(x,222,240)),
                                eccentricity = as.numeric(str_sub(x,184,202)),
                                meanAnomaly0 = as.numeric(str_sub(x,142,160)),
                                argumentOfPerigree  = as.numeric(str_sub(x,363,381)),
                                inclination0 = as.numeric(str_sub(x,325,343)),
                                omega0 = as.numeric(str_sub(x,283,301)),
                                meanMotionDifference  = as.numeric(str_sub(x,123,141)),
                                inclinationRate = as.numeric(str_sub(x,405,423)),
                                omegaRate = as.numeric(str_sub(x,382,400)),
                                cuc = as.numeric(str_sub(x,165,183)),
                                cus = as.numeric(str_sub(x,203,221)),
                                crc = as.numeric(str_sub(x,344,362)),
                                crs = as.numeric(str_sub(x,104,122)),
                                cic = as.numeric(str_sub(x,264,282)),
                                cis = as.numeric(str_sub(x,302,320)),
                                clockOffset = as.numeric(str_sub(x,24,42)),
                                clockDrift = as.numeric(str_sub(x,43,61)),
                                clockDriftRate = as.numeric(str_sub(x,62,80))),
             Cartesian = NA,
             # AODE = as.numeric(str_sub(x,85,103)),
             # # Codes = as.numeric(str_sub(x,424,442)),
             # BDTWeek = as.numeric(str_sub(x,443,461)),
             # # L2P_flag = as.numeric(str_sub(x,462,480)),
             # Accuracy = as.numeric(str_sub(x,485,503)),
             # SatH1 = as.numeric(str_sub(x,504,522)),
             # TGD1 = as.numeric(str_sub(x,523,541)),
             # TGD2 = as.numeric(str_sub(x,542,560)),
             # TransmissionTime = as.numeric(str_sub(x,565,583)),
             # AODC = as.numeric(str_sub(x,584,602))
        )
      )
    }
    
    return(switch(str_sub(x,1,1),
                  "C" = readBeidou(x),
                  "E" = readGalileo(x),
                  "G" = readGPS(x),
                  "R" = readGlonass(x),
                  stop("Unknown Constellation Type")
    ))
    
  }
  
  endOfHeader <- detect_index(rinex,~endsWith(.,"END OF HEADER"))
  data <- rinex[(endOfHeader+1):length(rinex)]
  messageStart <- which(!startsWith(data," "))
  messageEnd <- (messageStart[2:length(messageStart)]-1) %>% c(length(data))
  messages <- map2_chr(messageStart,messageEnd,~c(data[.x:.y]) %>% paste0(collapse=""))
  out <- map(messages,possibly(getMessage,NA))
  return(out)
}

KeplerToLocation <- function(time,kepler){
  U <- 3986004.418 *10^8 #WGS84
  omega_earth <- 7292115.0 * 10^-11 #WGS84
  SolveEccentric <- function(m,e){
    TOLERANCE <- 1e-6
    e1 <- m
    e2 <- m + e * sin(e1)
    while(abs(e2-e1)>TOLERANCE) {
      e1 <- e2
      e2 <- m + e * sin(e1)
    }
    return(e2)
  }
  
  time_k <- (time - kepler["toe"] + 302400) %% 604800 - 302400
  m_k <- m0 + (U^0.5 /(kepler["sqrtA"]^3) + kepler["deltaN"]) * time_k
  e_k <- SolveEccentric(m_k,kepler["e"])
  v_k <- atan((1-kepler["e"]^2)^0.5 * sin(e_k) /(cos(e_k)-kepler["e"]))
  theta_k <- kepler["w"]+v_k
  u_k <- kepler["omega"]+v_k + kepler["cuc"] * cos(2*theta_k) + kepler["cus"] * sin(2*theta_k)  
  r_k <- kepler["sqrtA"]^2 * (1-kepler["e"]* cos(e_k)) + kepler["crc"]*cos(2*theta_k) + kepler["crs"] * sin(2*theta_k)
  i_k <- kepler["i0"]+kepler["iDot"]*time_k + kepler["cic"]*cos(2*theta_k) + kepler["cis"] * sin(2*theta_k)  
  omega_k <- kepler["omega0"] + (kepler["omegaDot"] - omega_earth) * time_k - omega_earth * kepler["toe"]
  
  xdash_k <- r_k * cos(u_k)
  ydash_k <- r_k * sin(u_k)
  
  x_k <- xdash_k * cos(omega_k) - ydash_k * cos(i_k) * sin(omega_k)
  y_k <- xdash_k * sin(omega_k) + ydash_k * cos(i_k) * cos(omega_k)
  z_k <- ydash_k * sin(i_k)
  
  return(c(x=x_k,y=y_k,z=z_k))
}

as.Kepler <- function(timeOfEphemeris,sqrtMajorAxis,eccentricity,meanAnomaly0,
                      argumentOfPerigree,inclination0,omega0,meanMotionDifference,
                      inclinationRate,omegaRate,cuc,cus,crc,crs,cic,cis,
                      clockOffset,clockDrift,clockDriftRate){
  return(c(toe=timeOfEphemeris,sqrtA=sqrtMajorAxis,e=eccentricity,m0=meanAnomaly0,
           w=argumentOfPerigree,i0=inclination0,omega0=omega0,deltaN=meanMotionDifference,
           iDot=inclinationRate,omegaDot=omegaRate,cuc=cuc,cus=cus,crc=crc,crs=crs,cic=cic,cis=cis,
           clockOffset=clockOffset,clockDrift=clockDrift,clockDriftRate=clockDriftRate))
}

as.Cartesian <- function(timeOfEphemeris,x,y,z,v_x,v_y,v_z,a_x,a_y,a_z,clockOffset,clockFrequencyOffset){
  return(c(toe=timeOfEphemeris,x=x,y=y,z=z,v_x=v_x,v_y=v_y,v_z=v_z,a_x=a_x,a_y=a_y,a_z=a_z,clockOffset=clockOffset,clockFrequencyOffset=clockFrequencyOffset))
}

CartesianToLocation <- function(time,cart){
  x_a <- cart["x"] * cos(theta_g) - cart["y"] * sin(theta_g)
  y_a <- cart["x"] * csin(theta_g) + cart["y"] * cos(theta_g)
  z_a <- cart["z"]
  v_x_a <- cart["v_x"] * cos(theta_g) - cart["v_y"] * sin(theta_g) - omega_earth * y_a
  v_y_a <-cart["v_x"] * sin(theta_g) + cart["v_y"] * cos(theta_g) - omega_earth * x_a
  v_z_a <- cart["v_z"]
  
  return(NA)
}
