
# All scripts named 1AB refer to the same processing steps done in script 1A and 1B for studies downloaded directly from Movebank
# adapted for dataset that were sent externally (as csv) because of different data formats (different devices)

#____________________
# DAILY DIARIES ####
#____________________


library(data.table)
library(move)
library(sf)
library(amt)
library(lubridate)
library(geosphere)

setwd("...")

source("Scripts/COT_publ/COT_functions.R") #For direction360


#___________________________________
# Andean Condor (Emily Shepard) ####

# Import
options(digits=6)
options(digits.secs=4)

df <- read.csv("DataAvailable/GpsAcc_ExternalDatasets/EmilysData/Condors/allcondors_GPS_DD_weather_dir_marked.csv", as.is=T)

length(unique(df$condor))

# In the condor data, the vedba was averaged across 10 seconds (5 sec before and after each GPS point) at a frequency of 40 Hz
# Some of the Time values have a space at the beginning, remove it:
table(nchar(df$Time))
df$Time[nchar(df$Time)==9] <- substr(df$Time[nchar(df$Time)==9], 2, 9)
table(nchar(df$Time))
# Now we can recreate the column datetime and transform it in POSIXct
df$datetime <- paste0(as.Date(df$Date, format="%d/%m/%Y")," ",df$Time)
df$datetime <- as.POSIXct(df$datetime, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
df <- df[order(df$condor,df$datetime),]
# The column "flight type" can be ignored
df$flight.type <- NULL
# The column "marked.event" indicated the condor behaviour classified by Hannah: NA = non flight, any number > 1 = flight (2=gliding, 3=slope soaring, 4=thermal soaring, 5=flapping)
table(df$marked.event)
table(is.na(df$marked.event))
# Split by individual
df_ls <- split(df, df$condor)
# Check timelag
i=6
summary(as.numeric(difftime(df_ls[[i]]$datetime[-1],df_ls[[i]]$datetime[-nrow(df_ls[[i]])], units="secs")))

# Plot the tracks
cols <- rainbow(length(df_ls))
plot(df$Longitude, df$Latitude, type="n")
lapply(1:length(df_ls), function(i) lines(df_ls[[i]]$Longitude, df_ls[[i]]$Latitude, col=cols[i]))

# Calculate track geometry infos

gpsAcc_ls <- split(df, df$condor)

# table(table(gpsAcc$Longitude)>1)
# table(table(gpsAcc$Latitude)>1)
gpsAcc_ls <- lapply(gpsAcc_ls, function(gpsAcc){
  print(unique(gpsAcc$condor))
  # Order trip observations by timestamp
  gpsAcc <- gpsAcc[order(gpsAcc$datetime),]
  # Thin the data to 10+-5 min
  indAmt <- mk_track(tbl=gpsAcc, all_cols=T,
                     .x=Longitude, .y=Latitude, crs = st_crs(4326),
                     .t=datetime, order_by_ts = T, check_duplicates = T)
  indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
  m <- as_move(indAmt)
  # Remove outliers based on speed (max 50 m/s)
  while(any(move::speed(m) >= 50) == T){
    m <- m[which(c(NA, move::speed(m)) < 50)]
  }
  # Add variables about the track geometry
  m$timeLag_min <- c(NA, timeLag(m, units="mins"))
  m$altitudeDiff <- c(NA, (m$Altitude[-1] - m$Altitude[-nrow(m)])) 
  m$vertSpeed_ms <- m$altitudeDiff/(m$timeLag_min*60)
  m$stepLength_m <- c(NA, move::distance(m))
  m$groundSpeed_ms <- c(NA, move::speed(m))
  m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
  m$turnAngle <- c(NA, turnAngleGc(m), NA)  # Add variables about the track geometry
  m$individual.local.identifier <- m@idData$condor

  return(as.data.frame(m))
})

# Exclude empty elements from the list and bind it
gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
df_geom <- as.data.frame(rbindlist(gpsAcc_ls, fill=T))
# Import body mass infos and merge body mass info to the dataframe
speciesInfo <- read.csv("/home/mscacco/ownCloud/Martina/ProgettiVari/COT/DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)
df_geom$individual.taxon.canonical.name <- "Vultur gryphus"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
names(df_geom)[names(df_geom) %in% c("Altitude","timestamps")] <- c("height","timestamp")
names(df_geom)[names(df_geom) == "mean_VeDBA"] <- "meanVedba"
df_geom$study.name <- "Andean Condor Vultur gryphus Bariloche, Argentina, 2013-2018"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "DailyDiary"
# Add infos about ACC sampling
# In the condor data, the vedba was averaged across 10 seconds (5 sec before and after each GPS point) at a frequency of 40 Hz
df_geom$acc_sampl_freq_per_axis <- 40
df_geom$n_samples_per_axis <- 400
df_geom$acc_burst_duration_s <- 10
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba),]
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_EmilysDataHannah_Vultur-gryphus_DD_dfGpsAcc_geom_10min.RData")


#_____________________________________
# Tropicbirds (Emily Shepard) ####

# Ignore marked event
# Sampling frequency is 40 Hz
# VeDBA is averaged over a minute (30 s before and after each GPS point) that is on 60*40=2400 ACC observations

# Import
options(digits.secs=4)

df <- read.table("DataAvailable/GpsAcc_ExternalDatasets/EmilysData/TropicBirds/TropicBirds_Subset_GPS_FlightStyle_Coord2+VeDBA.txt", header=T, row.names = NULL, na.strings = c("NA",""), as.is=T)
df$timestamp <- paste(df$row.names, df$timestamp, sep=" ")
df$row.names <- NULL
df$timestamp <- as.POSIXct(df$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
#identical(df$TripID, df$TripID2)
df$TripID2 <- NULL
df$BirdID <- sapply(strsplit(df$TripID, "_"), "[", 1)
df <- df[order(df$BirdID,df$TripID,df$timestamp),]
#identical(df$Marked.event.x, df$Marked.event)
table(df$Marked.event.x)
table(df$Marked.event)
boxplot(VeDBA~Marked.event.x, data=df)

# This dataset contains only flight data:
table(df$Activity)
table(df$Label)

# Split by individual and trip
df_ls <- split(df, df$TripID)
# Check timelag
lapply(df_ls, function(df){
  summary(as.numeric(difftime(df$timestamp[-1],df$timestamp[-nrow(df)], units="secs")))
})

# Plot the tracks
cols <- rainbow(length(df_ls))
png("DataAvailable/GpsAcc_ExternalDatasets/EmilysData/TropicBirds/rawTraj_allBirds.png", width=8,height=8,units="in",res=300)
plot(df$location.long, df$location.lat, type="n")
lapply(1:length(df_ls), function(i) lines(df_ls[[i]]$location.long, df_ls[[i]]$location.lat, col=cols[i]))
dev.off()

# Calculate track geometry infos

gpsAcc_ls <- split(df, df$TripID)

gpsAcc_ls <- llply(gpsAcc_ls, function(gpsAcc){
  print(unique(gpsAcc$TripID))
  # Order trip observations by timestamp
  gpsAcc <- as.data.frame(gpsAcc[order(gpsAcc$timestamp),])
  # Thin the data to 10+-5 min
  indAmt <- mk_track(tbl=gpsAcc, all_cols=T,
                     .x=location.long, .y=location.lat, crs = st_crs(4326),
                     .t=timestamp, order_by_ts = T, check_duplicates = T)
  indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
  m <- as_move(indAmt)
  # Remove outliers based on speed (max 50 m/s)
  while(any(move::speed(m) >= 50) == T){
    m <- m[which(c(NA, move::speed(m)) < 50)]
  }
  # Add variables about the track geometry
  m$timeLag_min <- c(NA, timeLag(m, units="mins"))
  m$altitudeDiff <- c(NA, (m$Altitude[-1] - m$Altitude[-nrow(m)])) 
  m$vertSpeed_ms <- m$altitudeDiff/(m$timeLag_min*60)
  m$stepLength_m <- c(NA, move::distance(m))
  m$groundSpeed_ms <- c(NA, move::speed(m))
  m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
  m$turnAngle <- c(NA, turnAngleGc(m), NA)  # Add variables about the track geometry
  m$individual.local.identifier <- m@idData$BirdID
  m$track_flight_id <- m@idData$TripID

  return(as.data.frame(m))
}, .parallel=T)
# Exclude empty elements from the list
gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
# Import body mass infos
speciesInfo <- read.csv("/home/mscacco/ownCloud/Martina/ProgettiVari/COT/DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)
# Bind, add body mass value, change some column names to match other datasets and save
df_geom <- do.call(rbind, gpsAcc_ls)
df_geom$individual.taxon.canonical.name <- "Phaethon rubricauda"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
names(df_geom)[names(df_geom) %in% c("Altitude","timestamps")] <- c("height","timestamp")
names(df_geom)[names(df_geom) == "VeDBA"] <- "meanVedba"
df_geom$study.name <- "Red-tailed tropicbirds (Phaethon rubricauda) Round Island"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "DailyDiary"
# In the tropicbird data, the vedba was averaged across 1 minute (30 sec before and after each GPS point) at a frequency of 40 Hz
df_geom$acc_sampl_freq_per_axis <- 40
df_geom$n_samples_per_axis <- 2400
df_geom$acc_burst_duration_s <- 60
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba),]
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_EmilysDataBaptiste_Phaethon-rubricauda_DD_dfGpsAcc_geom_10min.RData")


#_________________________________________________________
# Atlantic yellow-nosed albatross (Stefan Schoombie) ####

# Import
library(readr)
options(digits=6)
options(digits.secs=4)

df <- read_csv("DataAvailable/GpsAcc_ExternalDatasets/StefanSchoombieData/ayna_all_out.zip")

length(unique(df$id))

df$datetime <- as.POSIXct(df$datetime, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
df <- df[order(df$id,df$datetime),]
df$gps_datetime <- as.POSIXct(df$gps_datetime, format="%Y-%m-%d %H:%M:%OS", tz="UTC")

# ACC data are collected continuously but GPS are not. 
# Check timelags between ACC observations, they are about 0-0.02 seconds
acc_ls <- split(df, df$id)
summary(unlist(sapply(acc_ls, function(df_ind){
  df_ind <- df_ind[order(df_ind$datetime),]
  return(as.numeric(difftime(df_ind$datetime[-1], df_ind$datetime[-length(df_ind$datetime)], units="secs")))
})))

# isolate only gps and remove duplicated gps_datetime, gps_lat, gps_long
# ACC Sampling frequency is 40 Hz
gpsDf <- df[complete.cases(df[,c("gps_lon","gps_lat")]),]
gpsDf <- gpsDf[order(gpsDf$gps_datetime),]
gpsDf <- gpsDf[!duplicated(gpsDf[,c("gps_datetime","gps_lon","gps_lat")]),]
# add empty ACC variables
gpsDf$meanVedba <- NA
gpsDf$burstDur <- NA
gpsDf$sampFreq <- NA
# split by individual
gps_ls <- split(gpsDf, gpsDf$id)
# Check timelags between GPS locations, they are about 1 h or more
summary(unlist(sapply(gps_ls, function(df_ind){
  df_ind <- df_ind[order(df_ind$gps_datetime),]
  return(as.numeric(difftime(df_ind$datetime[-1], df_ind$datetime[-length(df_ind$datetime)], units="mins")))
})))

# Plot the tracks
cols <- rainbow(length(gps_ls))
plot(gpsDf$gps_lon, gpsDf$gps_lat, type="n")
lapply(1:length(gps_ls), function(i){
  lines(gps_ls[[i]]$gps_lon, gps_ls[[i]]$gps_lat, col=cols[i])
  points(gps_ls[[i]]$gps_lon, gps_ls[[i]]$gps_lat, col=cols[i])})

### Thin data to 5 minutes, associate ACC 10 seconds around each GPS point and calculate track geometry infos
table(duplicated(gpsDf[,c("id","fls","gps_datetime")]))

names(gps_ls) %in% names(acc_ls)

gpsThin_ls <- lapply(gps_ls, function(gps){
  print(unique(gps$id))
  # Order trip observations by timestamp
  gps <- gps[order(gps$gps_datetime),]
  # remove duplicates
  gps <- gps[!duplicated(gps[,c("id","gps_datetime")]),]
  # Thin the data to 10+-5 min
  # (for eobs this step was already done in previous step before associating acc)
  indAmt <- mk_track(tbl=gps, all_cols=T,
                     .x=gps_lon, .y=gps_lat, crs = st_crs(4326),
                     .t=gps_datetime, order_by_ts = T, check_duplicates = T)
  indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
  m <- as_move(indAmt)
  # Remove outliers based on speed (max 50 m/s)
  while(any(move::speed(m) >= 50) == T){
    m <- m[which(c(NA, move::speed(m)) < 50)]
  }
  # Add variables about the track geometry
  if(n.locs(m)>1){
    m$timeLag_min <- c(NA, timeLag(m, units="mins"))
    m$height <- NA
    m$altitudeDiff <- NA 
    m$vertSpeed_ms <- NA
    m$stepLength_m <- c(NA, move::distance(m))
    m$groundSpeed_ms <- c(NA, move::speed(m))
    m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
    m$turnAngle <- c(NA, turnAngleGc(m), NA)  # Add variables about the track geometry
    return(as.data.frame(m))
  }
})
# Exclude empty elements and errors from the list and bind it
gpsThin_ls <- gpsThin_ls[which(!sapply(gpsThin_ls, is.null))]
is.error <- function(x) inherits(x, "try-error")
gpsThin_ls <- gpsThin_ls[!vapply(gpsThin_ls, is.error, logical(1))]

# For each individual average VeDBA over 10 sec (5 s before and after each GPS point) that is on (40Hz*5)*2=400 ACC observations
df_geom_ls <- lapply(gpsThin_ls, function(gps){
  # extract ACC data corresponding to the same individual
  dfID <- acc_ls[[unique(gps$id)]]
  
  for(i in 2:nrow(gps)){
    gpsTime <- gps[i,"timestamps"]
    accSub <- dfID[dfID$datetime >= gpsTime-5 & dfID$datetime <= gpsTime+5, c("datetime","Acc_x","Acc_y","Acc_z","VeDBA")] #subset of acc within 10 s around each gps point (+/- 5 sec)
    if(nrow(accSub > 0)){
      gps$burstDur[i] <- round(as.numeric(difftime(accSub$datetime[nrow(accSub)], accSub$datetime[1], units="secs")))
      gps$sampFreq[i] <- 40 #(original sampl freq of 40 Hz (like all DD tags)
      gps$meanVedba[i] <- mean(accSub$VeDBA, na.rm=T) # average vedba and associate to the gps point
      # accSub$vedba <- sqrt((accSub[,"Acc_x"]-mean(accSub[,"Acc_x"]))^2 + (accSub[,"Acc_y"]-mean(accSub[,"Acc_y"]))^2 + (accSub[,"Acc_z"]-mean(accSub[,"Acc_z"]))^2) #calculate vedba
      # gps$meanVedba[i] <- mean(accSub$vedba) # average vedba and associate to the gps point
    }
  }
  return(gps)
})
df_geom <- do.call(rbind, df_geom_ls)
# Import body mass infos and merge body mass info to the dataframe
speciesInfo <- read.csv("DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)

df_geom$individual.taxon.canonical.name <- "Thalassarche chlororhynchos"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$species_latin == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$species_latin == unique(df_geom$individual.taxon.canonical.name)]
df_geom$study.name <- "Atlantic yellow-nosed albatross (Thalassarche chlororhynchos) from Stefan Schoombie"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "DailyDiary"
# remove unnecessary columns
df_geom$datetime <- NULL
df_geom$g <- NULL
df_geom[,grep("Acc|Mag|VeSBA|VeDBA", names(df_geom), value=T)] <- NULL
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba),]
# change column names to match other datasets
names(df_geom)[names(df_geom) %in% c("id","timestamps")] <- c("timestamp","individual.local.identifier")
names(df_geom)[names(df_geom) %in% c("burstDur","sampFreq")] <- c("acc_burst_duration_s","acc_sampl_freq_per_axis")
df_geom$n_samples_per_axis <- df_geom$acc_sampl_freq_per_axis * df_geom$acc_burst_duration_s
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_EmilysDataSchoombie_Thalassarche-chlororhynchos_DD_dfGpsAcc_geom_10min.RData")


