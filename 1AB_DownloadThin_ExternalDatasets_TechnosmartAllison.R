
# All scripts named 1AB refer to the same processing steps done in script 1A and 1B for studies downloaded directly from Movebank
# adapted for dataset that were sent externally (as csv) because of different data formats (different devices)

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

# error function for try
is.error <- function(x) inherits(x, "try-error")


# Set as wd the path where you stored the content downloaded from the Edmond repository
setwd("...")

# Create folder to store processed data, study by study
dir.create("DataProcessed")

# store the path where you downloaded the scripts
codePath <- "..."

source(paste0(codePath, "COT_publ/COT_functions.R")) #For direction360


#____________________________________________
# Data Allison Patterson, pre-processed ####

# list the zip folders containing the raw data per species
zipFolds <- list.files("DataAvailable/GpsAcc_ExternalDatasets/AllisonPattersonsData/rawAccData/rawConvertedAccData", pattern=".zip", full.names=T)

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
  gpsDf$acc_axes <- "XYZ"
  gpsDf$acc_Naxes <- 3
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

