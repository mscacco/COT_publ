
# All scripts named 1AB refer to the same processing steps done in script 1A and 1B for studies downloaded directly from Movebank
# adapted for dataset that were sent externally (as csv) because of different data formats (different devices)

#_______________
# ORNITELA ####
#_______________


setwd("...")

library(data.table)
library(ggplot2)
library(dplyr)
library(plotly)
library(lubridate)
library(geosphere)
library(amt)
library(move)
library(sf)
library(plyr)
library(doParallel)
#detectCores()
doParallel::registerDoParallel(5)


options(digits.secs=3)

source("Scripts/COT_publ/COT_functions.R")

#_______________
# ORNITELA ####
#_______________

#_____________________________________
# Trumpeter Swans (David Wolfson) ####


#_____
# Work with GPS data

fls <- list.files("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan", pattern = "csv", full.names = T)

allGPS <- as.data.frame(rbindlist(lapply(fls, function(f){
  df <- read.csv(f)
  ind <- gsub(".*/|.csv","",f)
  if(length(unique(df$device_id))>1){stop(paste0("The file '",f,"' contains more than one device/individual"))}
  gps <- df[df$datatype == "GPS" & complete.cases(df[,c("Longitude","Latitude")]),1:15]
  gps$individual.local.identifier <- paste0(ind,"_",unique(gps$device_id))
  gps$UTC_datetime <- as.POSIXct(gps$UTC_datetime, "%Y-%m-%d %H:%M:%S", tz="UTC")
  gps <- gps[order(gps$UTC_datetime),]
  gps <- gps[!(gps$Latitude == 0 & gps$Longitude == 0),]
  gps$timeLag_s <- c(NA,as.numeric(difftime(gps$UTC_datetime[-1], gps$UTC_datetime[-nrow(gps)], units = "secs")))
  
  return(gps)
})))
rm(df)

# data are collected every 15 min
summary(allGPS$timeLag_s)/60
table(allGPS$individual.local.identifier)

# save all gps from all individuals together
saveRDS(allGPS, "DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/onlyGPS_allIndividuals_raw.rds")

#_____
# Work with ACC data (one individual at a time)

# Make sure plyr is not loaded before running the following, cause it masks group_by
detach("package:plyr", unload = TRUE)
dir.create("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/plots_accTest")

fls <- list.files("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan", pattern = "csv", full.names = T)

allACC_perBursts <- as.data.frame(rbindlist(lapply(fls, function(f){
  ind <- gsub(".*/|.csv","",f)
  df <- read.csv(f)
  acc <- df[df$datatype == "SENSORS",c("device_id","UTC_datetime","acc_x"       ,"acc_y","acc_z")]
  acc$individual.local.identifier <- paste0(ind,"_",unique(acc$device_id))
  acc$UTC_datetime <- as.POSIXct(acc$UTC_datetime, "%Y-%m-%d %H:%M:%OS", tz="UTC")
  # inspect burst duration and sampl. freq
  timeLag_s <- c(0,as.numeric(difftime(acc$UTC_datetime[-1], acc$UTC_datetime[-nrow(acc)], units = "secs")))

  # ACC is collected in continuous bursts of 30 samples each, with a duration of 2.5 seconds, every 5 min (300 s) [~ 12 Hz]
  acc$burstID <- cumsum(timeLag_s > 1)

  # Keep regular bursts, 30 samples 12 Hz
  burstsToKeep <- names(table(acc$burstID)[table(acc$burstID) == 30])
  acc <- acc[acc$burstID %in% burstsToKeep,]
  
  # transform ACC axes in Gs
  acc$acc_xG <- acc$acc_x/1000
  acc$acc_yG <- acc$acc_y/1000
  acc$acc_zG <- acc$acc_z/1000
  # calculate vedba per sample (each ACC entry)
  acc$vedba <- sqrt((acc[,"acc_xG"]-mean(acc[,"acc_xG"]))^2 + (acc[,"acc_yG"]-mean(acc[,"acc_yG"]))^2 + (acc[,"acc_zG"]-mean(acc[,"acc_zG"]))^2) 
  
  # example plots
  set.seed(100)
  randomBursts <- sample(unique(acc$burstID),20)
  accSub <- acc[acc$burstID %in% randomBursts,]
  acc_plotLs <- split(accSub, accSub$burstID)

  lapply(acc_plotLs, function(x){
    # Ornitela doesn't show milliseconds so we can't order/plot by timestamps.
    # We need to trust and plot based on order in which samples occur in the burst
    png(paste0("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/plots_accTest/device_",unique(x$individual.local.identifier),"_burst_",unique(x$burstID),".png"))
    plot(1:nrow(x), x$acc_xG, type="n", ylim=c(-2.1,2.1))
    lines(1:nrow(x), x$acc_xG, col="red")
    lines(1:nrow(x), x$acc_yG, col="green")
    lines(1:nrow(x), x$acc_zG, col="blue")
    dev.off()
  })
  
  # summarise info per burst id (each is 2.5 s long)
  accPerBurst <- acc %>% group_by(burstID) %>% 
    summarise(individual.local.identifier=unique(individual.local.identifier),
              UTC_datetime_start=first(UTC_datetime),
              UTC_datetime_end=last(UTC_datetime),
              meanVedba_Gs=mean(vedba, na.rm=T),
              cumVedba_Gs=sum(vedba, na.rm=T))

  # add ACC burst sampling info found out above
  accPerBurst$acc_burst_duration_s <- 2.5
  accPerBurst$acc_sampl_freq_per_axis <- 12
  accPerBurst$n_samples_per_axis <- 30
  
  return(accPerBurst)
})))

hist(allACC_perBursts$meanVedba_Gs)
summary(allACC_perBursts$meanVedba_Gs)
length(unique(allACC_perBursts$individual.local.identifier)) #how many tags
table(allACC_perBursts$individual.local.identifier) #how many ACC bursts per tag

saveRDS(allACC_perBursts, "DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/onlyACC_allIndividuals_rawPerBurst.rds")

#____
# Associate ACC to GPS

allGPS <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/onlyGPS_allIndividuals_raw.rds")
allACC <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/onlyACC_allIndividuals_rawPerBurst.rds")

allGPS_ls <- split(allGPS, allGPS$individual.local.identifier)

accCols <- c("meanVedba_Gs", "cumVedba_Gs", "acc_burst_duration_s", "acc_sampl_freq_per_axis", "n_samples_per_axis")
timeTol <- 5*60 #in seconds

all_GpsAcc <- as.data.frame(rbindlist(lapply(allGPS_ls, function(gps){
  acc <- allACC[allACC$individual.local.identifier == unique(gps$individual.local.identifier),]
  if(nrow(acc)>0){
    gpsAcc <- ACCtoGPS(acc, gps,
                       timeTolerance=timeTol,
                       ColsToAssociate=accCols,
                       ACCtimeCol="UTC_datetime_start", GPStimeCol="UTC_datetime")
    return(gpsAcc)
  }
})))

table(is.na(all_GpsAcc$acc_event_id)) #how many gps points didn't get an ACC associated?
summary(all_GpsAcc$diff_acc_time_s) # what is the average time diff between gps and associated acc?

saveRDS(all_GpsAcc, "DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/GPS-ACC_allIndividuals_rawPerBurst.rds")

#______
# Calculate track geometry infos

all_GpsAcc <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/DavidWolfson_TrumpeterSwan/GPS-ACC_allIndividuals_rawPerBurst.rds")

gpsAcc_ls <- split(all_GpsAcc, all_GpsAcc$individual.local.identifier)

gpsAcc_ls <- llply(gpsAcc_ls, function(gpsAcc){
  print(unique(gpsAcc$individual.local.identifier))
  # Order individual observations by timestamp
  gpsAcc <- as.data.frame(gpsAcc[order(gpsAcc$UTC_datetime),])
  # Thin the data to 10+-5 min
  indAmt <- mk_track(tbl=gpsAcc, all_cols=T,
                     .x=Longitude, .y=Latitude, crs = st_crs(4326),
                     .t=UTC_datetime, order_by_ts = T, check_duplicates = T)
  indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
  m <- as_move(indAmt)
  # Remove outliers based on speed (max 50 m/s)
  while(any(move::speed(m) >= 50) == T){
    m <- m[which(c(NA, move::speed(m)) < 50)]
  }
  # Add variables about the track geometry
  m$timeLag_min <- c(NA, timeLag(m, units="mins"))
  m$altitudeDiff <- c(NA, (m$Altitude_m[-1] - m$Altitude_m[-nrow(m)])) 
  m$vertSpeed_ms <- m$altitudeDiff/(m$timeLag_min*60)
  m$stepLength_m <- c(NA, move::distance(m))
  m$groundSpeed_ms <- c(NA, move::speed(m))
  m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
  m$turnAngle <- c(NA, turnAngleGc(m), NA)  # Add variables about the track geometry
  m$tag.local.identifier <- m@idData$device_id
  
  return(as.data.frame(m))
}, .parallel=T)
# Exclude empty elements from the list
gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
# Import body mass infos
speciesInfo <- read.csv("DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)
# Bind, add body mass value, change some column names to match other datasets and save
df_geom <- rbindlist(gpsAcc_ls)
df_geom$individual.taxon.canonical.name <- "Cygnus buccinator"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
names(df_geom)[names(df_geom) %in% c("Altitude_m","timestamps")] <- c("height","timestamp")
df_geom$study.name <- "Trumpeter swan (Cygnus buccinator)"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "Ornitela"
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba_Gs),]
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_DavidWolfson_Cygnus-buccinator_Ornitela_dfGpsAcc_geom_10min.RData")


#______________________________
# Ospreys (Flavio Monti) ####

#_____
# Work with GPS data (all individuals)
fls <- list.files("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela", pattern = "csv", full.names = T)

allGPS <- as.data.frame(rbindlist(lapply(fls, function(f){
  ind <- gsub(".*ring|.csv","",f)
  df <- read.csv(f)
  if(length(unique(df$device_id))>1){stop(paste0("The file '",f,"' contains more than one device/individual"))}
  gps <- df[df$datatype == "GPS" & complete.cases(df[,c("Longitude","Latitude")]),1:15]
  gps$individual.local.identifier <- paste0(ind,"_",unique(gps$device_id))
  gps$UTC_datetime <- as.POSIXct(gps$UTC_datetime, "%Y-%m-%d %H:%M:%S", tz="UTC")
  gps <- gps[order(gps$UTC_datetime),]
  gps <- gps[!(gps$Latitude == 0 & gps$Longitude == 0),]
  gps$timeLag_s <- c(NA,as.numeric(difftime(gps$UTC_datetime[-1], gps$UTC_datetime[-nrow(gps)], units = "secs")))
  
  return(gps)
})))
rm(df)

# data are collected every 5 min, followed by a 20 sec bursts at 1 Hz
summary(allGPS$timeLag_s)/60
allGPS$timeLag_s[1000:1200]
length(unique(allGPS$individual.local.identifier))
table(allGPS$individual.local.identifier)

saveRDS(allGPS, "DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/onlyGPS_allIndividuals_raw.rds")

#_____
# Work with ACC data (one individual at a time)

# Make sure plyr is not loaded before running the following, as it masks group_by from dplyr
detach("package:plyr", unload = TRUE)

dir.create("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/plots_accTest")

fls <- list.files("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela", pattern = "csv", full.names = T)

allACC_perBursts <- as.data.frame(rbindlist(lapply(fls, function(f){
  df <- read.csv(f)
  ind <- gsub(".*ring|.csv","",f)
  acc <- df[df$datatype == "SENSORS",c("device_id","UTC_datetime","acc_x","acc_y","acc_z")]
  acc$individual.local.identifier <- paste0(ind,"_",unique(acc$device_id))
  acc$UTC_datetime <- as.POSIXct(acc$UTC_datetime, "%Y-%m-%d %H:%M:%OS", tz="UTC")
  # inspect burst duration and sampl. freq
  timeLag_s <- c(0,as.numeric(difftime(acc$UTC_datetime[-1], acc$UTC_datetime[-nrow(acc)], units = "secs")))
  # head(acc, 50)
  # summary(timeLag_s[timeLag_s>1]/60)
  # table(timeLag_s)
  
  # ACC is collected every 14 min (891 s) in continuous bursts of either 100 or 200 or 400 or 1200 samples each, with a duration of either 5 or 20 or 60 seconds [meaning always ~ 20 Hz]
  acc$burstID <- cumsum(timeLag_s > 1)
  # table(acc$burstID)
  # table(table(acc$burstID))
  # length(unique(acc$burstID))
  # ls <- split(acc, acc$burstID)
  # table(sapply(ls, function(b){as.numeric(difftime(max(b$UTC_datetime),min(b$UTC_datetime),units="secs"))}))

  # Keep regular bursts, the 100, 200 and 400 samples ones are by far the most common, so we keep those
  burstsToKeep <- names(table(acc$burstID)[table(acc$burstID) %in% c(100,200,400)])
  acc <- acc[acc$burstID %in% burstsToKeep,]
  
  # transform ACC axes in Gs
  acc$acc_xG <- acc$acc_x/1000
  acc$acc_yG <- acc$acc_y/1000
  acc$acc_zG <- acc$acc_z/1000
  # calculate vedba per sample (each ACC entry)
  acc$vedba <- sqrt((acc[,"acc_xG"]-mean(acc[,"acc_xG"]))^2 + (acc[,"acc_yG"]-mean(acc[,"acc_yG"]))^2 + (acc[,"acc_zG"]-mean(acc[,"acc_zG"]))^2) 
  
  # example plots
  set.seed(100)
  randomBursts <- sample(unique(acc$burstID),20)
  accSub <- acc[acc$burstID %in% randomBursts,]
  acc_plotLs <- split(accSub, accSub$burstID)
  
  lapply(acc_plotLs, function(x){
    # Ornitela doesn't show milliseconds so we can't order/plot by timestamps.
    # We need to trust and plot based on order in which samples occur in the burst
    png(paste0("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/plots_accTest/device_",unique(x$individual.local.identifier),"_burst_",unique(x$burstID),".png"))
    plot(1:nrow(x), x$acc_xG, type="n", ylim=c(-2.1,2.1))
    lines(1:nrow(x), x$acc_xG, col="red")
    lines(1:nrow(x), x$acc_yG, col="green")
    lines(1:nrow(x), x$acc_zG, col="blue")
    dev.off()
  })
  
  # summarise info per burst id (they have different duration)
  accPerBurst <- acc %>% group_by(burstID) %>% 
    summarise(individual.local.identifier=unique(individual.local.identifier),
              UTC_datetime_start=first(UTC_datetime),
              UTC_datetime_end=last(UTC_datetime),
              meanVedba_Gs=mean(vedba, na.rm=T),
              cumVedba_Gs=sum(vedba, na.rm=T),
              n_samples_per_axis=length(UTC_datetime))
  
  # add ACC burst sampling info found out above
  accPerBurst$acc_burst_duration_s <- as.numeric(difftime(accPerBurst$UTC_datetime_end, accPerBurst$UTC_datetime_start, units="secs"))
  accPerBurst$acc_burst_duration_s[accPerBurst$acc_burst_duration_s==21] <- 20
  accPerBurst$acc_sampl_freq_per_axis <- round(accPerBurst$n_samples_per_axis/accPerBurst$acc_burst_duration_s)
  
  return(accPerBurst)
})))

head(allACC_perBursts)
table(allACC_perBursts$acc_sampl_freq_per_axis)
allACC_perBursts <- allACC_perBursts[allACC_perBursts$acc_sampl_freq_per_axis %in% c(20,25),]
table(allACC_perBursts$acc_sampl_freq_per_axis)
table(allACC_perBursts$acc_burst_duration_s)
table(allACC_perBursts$n_samples_per_axis)
hist(allACC_perBursts$meanVedba_Gs)
summary(allACC_perBursts$meanVedba_Gs)
length(unique(allACC_perBursts$individual.local.identifier)) #how many individual-tag combinations
table(allACC_perBursts$individual.local.identifier) #how many ACC bursts per ind

saveRDS(allACC_perBursts, "DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/onlyACC_allIndividuals_rawPerBurst.rds")

#___
# Associate ACC to GPS

allGPS <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/onlyGPS_allIndividuals_raw.rds")
allACC <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/onlyACC_allIndividuals_rawPerBurst.rds")

allGPS_ls <- split(allGPS, allGPS$individual.local.identifier)

accCols <- c("meanVedba_Gs", "cumVedba_Gs", "acc_burst_duration_s", "acc_sampl_freq_per_axis", "n_samples_per_axis")
timeTol <- 5*60 #in seconds

all_GpsAcc <- as.data.frame(rbindlist(lapply(allGPS_ls, function(gps){
  acc <- allACC[allACC$individual.local.identifier == unique(gps$individual.local.identifier),]
  if(nrow(acc)>0){
    gpsAcc <- ACCtoGPS(acc, gps,
                       timeTolerance=timeTol,
                       ColsToAssociate=accCols,
                       ACCtimeCol="UTC_datetime_start", GPStimeCol="UTC_datetime")
    return(gpsAcc)
  }
})))

table(is.na(all_GpsAcc$acc_event_id)) #how many gps points didn't get an ACC associated?
summary(all_GpsAcc$diff_acc_time_s) # what is the average time diff between gps and associated acc?
length(unique(all_GpsAcc$individual.local.identifier))

saveRDS(all_GpsAcc, "DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/GPS-ACC_allIndividuals_rawPerBurst.rds")


#___
# Calculate track geometry infos

all_GpsAcc <- readRDS("DataAvailable/GpsAcc_ExternalDatasets/FlavioMonti_osprey_ornitela/GPS-ACC_allIndividuals_rawPerBurst.rds")

# Split dataset by individual/device
gpsAcc_ls <- split(all_GpsAcc, all_GpsAcc$individual.local.identifier)

gpsAcc_ls <- lapply(gpsAcc_ls, function(gpsAcc)try({
  print(unique(gpsAcc$individual.local.identifier))
  # Order individual observations by timestamp
  gpsAcc <- as.data.frame(gpsAcc[order(gpsAcc$UTC_datetime),])
  # Remove duplicated times (they occur in high resolution GPS bursts)
  if(length(which(duplicated(gpsAcc$UTC_datetime)))>0){
    gpsAcc <- gpsAcc[-which(duplicated(gpsAcc$UTC_datetime)),]}
  # Thin the data to 10+-5 min
  indAmt <- mk_track(tbl=gpsAcc, all_cols=T,
                     .x=Longitude, .y=Latitude, crs = st_crs(4326),
                     .t=UTC_datetime, order_by_ts = T, check_duplicates = T)
  indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
  m <- as_move(indAmt)
  # Remove outliers based on speed (max 50 m/s)
  while(any(move::speed(m) >= 50) == T){
    m <- m[which(c(NA, move::speed(m)) < 50)]
  }
  # Add variables about the track geometry
  m$timeLag_min <- c(NA, timeLag(m, units="mins"))
  m$altitudeDiff <- c(NA, (m$Altitude_m[-1] - m$Altitude_m[-nrow(m)])) 
  m$vertSpeed_ms <- m$altitudeDiff/(m$timeLag_min*60)
  m$stepLength_m <- c(NA, move::distance(m))
  m$groundSpeed_ms <- c(NA, move::speed(m))
  m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
  m$turnAngle <- c(NA, turnAngleGc(m), NA)  # Add variables about the track geometry
  m$tag.local.identifier <- m@idData$device_id
  
  return(as.data.frame(m))
}))#, .parallel=T)
# Exclude empty elements from the list
gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
# Import body mass infos
speciesInfo <- read.csv("DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)
# Bind, add body mass value, change some column names to match other datasets and save
df_geom <- rbindlist(gpsAcc_ls)
df_geom$individual.taxon.canonical.name <- "Pandion haliaetus"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$matchingSpeciesName == unique(df_geom$individual.taxon.canonical.name)]
names(df_geom)[names(df_geom) %in% c("Altitude_m","timestamps")] <- c("height","timestamp")
df_geom$study.name <- "Ospreys (Pandion haliaetus)"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "Ornitela"
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba_Gs),]
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_FlavioMonti_Pandion-haliaetus_Ornitela_dfGpsAcc_geom_10min.RData")

