#____________________
# TECHNOSMART ####
#____________________


library(data.table)
library(lubridate)
library(amt)
library(move)
library(geosphere)
library(ggplot2)
theme_set(theme_bw())

library(sf)
library(maps)
world_map <- map_data("world")

options(digits.secs=4)

library(plyr)
library(doParallel)
doParallel::registerDoParallel(6)

setwd("...")
source("Scripts/COT_publ/COT_functions.R")

# error function for try
is.error <- function(x) inherits(x, "try-error")


#____________________________________________
# Data Allison Patterson, pre-processed ####

# list the zip folders containing the raw data per species
zipFolds <- list.files("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData/rawConvertedAccData_csv_5_species", pattern=".zip", full.names=T)
# Import the deployment lists with species name and start/end deployments 
dep_data_sub <- read.csv("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/pre-processedData/AllisonsData_depId-species4accIndividuals_removedExclude.csv", as.is=T)
dep_data_sub$time_released <- as.POSIXct(dep_data_sub$time_released, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
dep_data_sub$time_recaptured <- as.POSIXct(dep_data_sub$time_recaptured, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
anyNA(dep_data_sub$time_released)
anyNA(dep_data_sub$time_recaptured)
# From the deployments table create a dataset with matching species code and scientific name
SpeciesIDs <- aggregate(species_sci~species, data=dep_data_sub, FUN=unique)

# Import body mass infos and check that all names match
speciesInfo <- read.csv("/DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)
table(unique(SpeciesIDs$species_sci) %in% speciesInfo$matchingSpeciesName)

# create directory for plots
dir.create("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData/TrajectoryPlots")

# Now run the following for each of the 5 species and each file, to associate vedba values to gps
for(z in 1:nrow(SpeciesIDs)){
  spCode <- SpeciesIDs$species[z]
  spName <- SpeciesIDs$species_sci[z]
  # subset the deployments
  deplSub <- dep_data_sub[which(dep_data_sub$species == spCode),]
  # unzip the content of the folder in a temporary file and list its content
  tmpdir <- paste(tempdir(), spCode, sep="/")
  unzip(zipFolds[grep(spCode, zipFolds)], exdir=tmpdir)
  fls <- list.files(tmpdir, pattern="\\.csv", full.names = T)
  # Of those files, keep only suitable deployments (dep_id in the deplSub dataframe)
  filesID <- gsub("A.[0-9]_|_A.[0-9]|A.[0-9] ", "", gsub(".*/|_S.*", "", fls))
  deplSub$dep_id_match <- gsub("A.[0-9]_|_A.[0-9]|A.[0-9] ", "", gsub("_S.*", "", deplSub$dep_id))
  fls_sub <- fls[filesID %in% deplSub$dep_id_match]
  filesID <- gsub("A.[0-9]_|_A.[0-9]|A.[0-9] ", "", gsub(".*/|_S.*", "", fls_sub))
  # For each file, separate GPS from ACC and associate to each GPS point the VeDBA calculated 5 sec before and after each GPS point
  gpsAcc_ls <- llply(1:length(filesID), function(j)try({
    print(filesID[j])
    f <- fls_sub[j]
    # For file n.7 of RBME we need to remove some rows, which are a replicate:
    if(spCode=="RBME" & filesID[j]=="13RBME20210723"){
      tmpFl <- readLines(f)
      tmpFl <- prov[-c(14136412:(14136412+49))]
      df <- read.table(text = tmpFl, header=T, sep="\t", na.strings = c("NA",""), as.is=T) #since readLines is a vector we need the argument text=
      df$Timestamp <- as.POSIXct(df$Timestamp, format = "%Y/%m/%d %H:%M:%OS", tz="UTC")
    }else if(spCode=="TBMU" & filesID[j] %in% c("A1 118600205 20190730","99683077 20190712","117639829 20190707")){
      df <- read.table(f, header=T, sep=",", na.strings = c("NA",""), as.is=T)
      df$Timestamp <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%OS", tz="UTC")
    }else{
      df <- read.table(f, header=T, sep="\t", na.strings = c("NA",""), as.is=T)
      df$Timestamp <- as.POSIXct(df$Timestamp, format = "%Y/%m/%d %H:%M:%OS", tz="UTC")
    }
    if(anyNA(df$Timestamp)){warning("There are NAs in the timestamp column.")}
    df <- df[order(df$TagID,df$Timestamp),]
    # Add the dep_id to the df for matching it with the deployment dataset
    df$dep_id_match <- filesID[j]
    # Subset the df based on the deployment start and end
    deplTag <- deplSub[deplSub$dep_id_match == unique(df$dep_id_match),]
    if(is.na(deplTag$time_recaptured)){
      df <- df[df$Timestamp > deplTag$time_released,]
    }else{df <- df[df$Timestamp > deplTag$time_released & df$Timestamp < deplTag$time_recaptured,]}
    # Extract only the gps point
    gps <- df[complete.cases(df[,c("location.lat","location.lon")]),]
    gps <- gps[order(gps$Timestamp),]
    # associate the mean VeDBA averaged 10 sec around each gps point
    gps$meanVedba <- NA
    gps$burstDur <- NA
    gps$sampFreq <- NA
    for(i in 1:nrow(gps)){
      gpsTime <- gps[i,"Timestamp"]
      accSub <- df[df$Timestamp > gpsTime-5 & df$Timestamp < gpsTime+5, c("Timestamp","X","Y","Z")] #subset of acc within 10 s around each gps point (+/- 5 sec)
      gps$burstDur[i] <- round(as.numeric(difftime(accSub$Timestamp[nrow(accSub)], accSub$Timestamp[1])))
      gps$sampFreq[i] <- round(nrow(accSub)/gps$burstDur[i], 1) #(sampl freq of 25 or 50 Hz)
      accSub$vedba <- sqrt((accSub[,"X"]-mean(accSub[,"X"]))^2 + (accSub[,"Y"]-mean(accSub[,"Y"]))^2 + (accSub[,"Z"]-mean(accSub[,"Z"]))^2) #calculate vedba
      gps$meanVedba[i] <- mean(accSub$vedba) # average vedba and associate to the gps point
    }
    return(gps)
  }), .parallel=T)
  # Remove potential individuals that returned errors during download
  gpsAcc_ls <- gpsAcc_ls[!vapply(gpsAcc_ls, is.error, logical(1))]
  # Exclude empty elements from the list
  gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
  # Bind all gps data with assocaited vedba, from all individuals, in one object
  gpsDf <- as.data.frame(rbindlist(gpsAcc_ls, fill=T))
  # Plot trajectories
  cols <- rainbow(length(gpsAcc_ls))
  png(paste0("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData/TrajectoryPlots/rawTraj_allBirds_",spCode,".png"), width=8,height=8,units="in",res=300)
  plot(location.lat~location.long, data=gpsDf, type="n", main=paste(spCode, spName, sep=" - "))
  lapply(1:length(gpsAcc_ls), function(x){lines(location.lat~location.long, data=gpsAcc_ls[[x]], col=cols[[x]])})
  dev.off()
  # add some column for later on binding this species with the other
  gpsDf$individual.taxon.canonical.name <- spName
  gpsDf$species_english <- speciesInfo$species_english[speciesInfo$matchingSpeciesName == spName]
  gpsDf$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$matchingSpeciesName == spName]
  names(gpsDf)[names(gpsDf) %in% c("TagID","Timestamp","location.lon")] <- c("individual.local.identifier","timestamp","location.long")
  gpsDf$study.name <- if(spCode=="RBME"){paste0("KyleElliott_",spCode)}else{paste0("AllisonPatterson_",spCode)}
  gpsDf$event.id <- 1:nrow(gpsDf)
  gpsDf$deviceType <- "Technosmart_axytreck"
  # In this data, we averaged vedba across 10 seconds (5 sec before and after each GPS point) at a frequency of 25 or 50 Hz
  names(gpsDf)[names(gpsDf) %in% c("burstDur","sampFreq")] <- c("acc_burst_duration_s","acc_sampl_freq_per_axis")
  gpsDf$n_samples_per_axis <- gpsDf$acc_sampl_freq_per_axis * gpsDf$acc_burst_duration_s
  # Finally save the dataset for this species
  save(gpsDf, file=paste0("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData/studyId_AllisonsData_",spCode,"_",sub(" ","-",spName),"_Technosmart_gpsAcc.RData"))
}

# Calculate track geometry

fls <- list.files("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData", pattern="AllisonsData.*.RData", full.names = T)

lapply(fls, function(f){
  load(f) #object gpsDf
  spName <- unique(gpsDf$individual.taxon.canonical.name)
  # Split dataset by tag id
  df_ls <- split(gpsDf, gpsDf$individual.local.identifier)
  
  # Track geometry
  gpsAcc_ls <- lapply(df_ls, function(gpsAcc){
    print(unique(gpsAcc$individual.local.identifier))
    # Order trip observations by timestamp
    gpsAcc <- as.data.frame(gpsAcc[order(gpsAcc$timestamp),])
    # Thin the data to 5 min
    #gpsAcc <- gpsAcc[!duplicated(gpsAcc[,c("timestamp","location.long","location.lat")]),]
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
    return(as.data.frame(m))
  })
  # Exclude empty elements from the list
  gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
  # Bind and save dataset for next step
  df_geom <- do.call(rbind, gpsAcc_ls)
  save(df_geom, file=paste0("DataProcessed/studyId_AllisonsData_",sub(" ","-",spName),"_Technosmart_dfGpsAcc_geom_10min.RData"))
  
})


#___________________________________________
# Barn own (Paolo Becciu) pre-processed ####

df <- read.csv("DataAvailable/GpsAcc_ExternalDatasets/PaoloBecciuData/barnOwls_annotated_individuals_2020_complete.csv", as.is=T)

# ACC data are collected continuously at 50 Hz but are then averaged per second to match the GPS sampling frequency. 

df$timestamp <- as.POSIXct(df$timestamp, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
# remove missing coorinates
df <- df[complete.cases(df[,c("location.lon","location.lat")]),]
# add empty ACC variables
df$meanVedba <- NA
df$burstDur <- NA
df$sampFreq <- NA
# split by individual
gps_ls <- split(df, df$TagID)
# Check timelags between GPS observations, they are about 1 second
summary(unlist(sapply(gps_ls, function(df_ind){
  df_ind <- df_ind[order(df_ind$timestamp),]
  return(as.numeric(difftime(df_ind$timestamp[-1], df_ind$timestamp[-length(df_ind$timestamp)], units="secs")))
})))
# Plot the tracks
cols <- rainbow(length(gps_ls))
plot(df$location.lon, df$location.lat, type="n")
lapply(1:length(gps_ls), function(i) lines(gps_ls[[i]]$location.lon, gps_ls[[i]]$location.lat, col=cols[i]))

# Thin data to 5 minutes, associate ACC 10 seconds around each GPS point and calculate track geometry infos

table(duplicated(df[,c("TagID","timestamp")]))

gpsThin_ls <- lapply(gps_ls, function(gpsAcc)try({
  print(unique(gpsAcc$TagID))
  # Order trip observations by timestamp
  gpsAcc <- gpsAcc[order(gpsAcc$timestamp),]
  # remove duplicates
  gpsAcc <- gpsAcc[!duplicated(gpsAcc[,c("TagID","timestamp")]),]
  # Thin the data to 10+-5 min
  # (for eobs this step was already done in previous step before associating acc)
  indAmt <- mk_track(tbl=gpsAcc, all_cols=T,
                     .x=location.lon, .y=location.lat, crs = st_crs(4326),
                     .t=timestamp, order_by_ts = T, check_duplicates = T)
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
}))
# Exclude empty elements and errors from the list and bind it
gpsThin_ls <- gpsThin_ls[which(!sapply(gpsThin_ls, is.null))]
is.error <- function(x) inherits(x, "try-error")
gpsThin_ls <- gpsThin_ls[!vapply(gpsThin_ls, is.error, logical(1))]
# For each individual average VeDBA over 10 sec (5 s before and after each GPS point) that is on 60*1=60 ACC observations
df_geom_ls <- lapply(gpsThin_ls, function(gps){
  dfID <- df[df$TagID == unique(gps$TagID),]
  for(i in 2:nrow(gps)){
    gpsTime <- gps[i,"timestamps"]
    accSub <- dfID[dfID$timestamp >= gpsTime-5 & dfID$timestamp <= gpsTime+5, c("timestamp","X","Y","Z","VeDBA")] #subset of acc within 10 s around each gps point (+/- 5 sec)
    if(nrow(accSub > 0)){
      gps$burstDur[i] <- round(as.numeric(difftime(accSub$timestamp[nrow(accSub)], accSub$timestamp[1], units="secs")))
      gps$sampFreq[i] <- 50 #(original sampl freq of 50 Hz, reduced to a 1 Hz VeDBA average)
      gps$meanVedba[i] <- mean(accSub$VeDBA) # average vedba and associate to the gps point
    }
  }
  return(gps)
})
df_geom <- do.call(rbind, df_geom_ls)
# Import body mass infos and merge body mass info to the dataframe
speciesInfo <- read.csv("/home/mscacco/ownCloud/Martina/ProgettiVari/COT/DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", as.is=T)

df_geom$individual.taxon.canonical.name <- "Tyto alba"
df_geom$species_english <- speciesInfo$species_english[speciesInfo$species_latin == unique(df_geom$individual.taxon.canonical.name)]
df_geom$BodyMass_value <- speciesInfo$BodyMass_value[speciesInfo$species_latin == unique(df_geom$individual.taxon.canonical.name)]
df_geom$study.name <- "Barn Ownls (Tyto alba) from P.Becciu and K.Schalcher"
df_geom$event.id <- 1:nrow(df_geom)
df_geom$deviceType <- "Technosmart"
# remove unnecessary columns
df_geom[,grep("dynamic|static|ODBA|VeDBA|prey", names(df_geom), value=T)] <- NULL
# change column names to match other datasets
names(df_geom)[names(df_geom) %in% c("TagID","timestamps")] <- c("timestamp","individual.local.identifier")
names(df_geom)[names(df_geom) %in% c("burstDur","sampFreq")] <- c("acc_burst_duration_s","acc_sampl_freq_per_axis")
df_geom$n_samples_per_axis <- df_geom$acc_sampl_freq_per_axis * df_geom$acc_burst_duration_s
# Remove data points for which no ACC information is available
df_geom <- df_geom[complete.cases(df_geom$meanVedba),]
# save the updated dataset with geometry and additional columns
save(df_geom, file="DataProcessed/studyId_PaoloBecciuData_Tyto-alba_Technosmart_dfGpsAcc_geom_10min.RData")






