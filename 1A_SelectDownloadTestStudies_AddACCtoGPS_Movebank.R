#____________________________________________________
# Process datasets/studies available on Movebank ####
#____________________________________________________

setwd("...")

#____________________________________________________________
# Import general info from studies identified in step 0A ####

allStudies <- read.csv("DataAvailable/GpsAcc_MovebankDatasets/allACCstudies_BirdsBats_downloadPermissionT_old.csv", as.is=T)

allStudies$timestamp_first_deployed_location <- strptime(allStudies$timestamp_first_deployed_location, format="%Y-%m-%d %H:%M:%S", tz="UTC")
allStudies$timestamp_last_deployed_location <- strptime(allStudies$timestamp_last_deployed_location, format="%Y-%m-%d %H:%M:%S", tz="UTC")
allStudies[,c("species_english","species","number_of_deployed_locations","deployment_duration_days","study_name","BodyMass_value")]

# Filter out studies with very few locations
# Let's say that we want a minimum of 1 week worth of data collected at 15 min sampling
minLocs <- 7*24*4 # 4 obs per hour, per 7 days
studiesSub <- allStudies[which(allStudies$number_of_deployed_locations > minLocs &
                           allStudies$deployment_duration_days > 7),]
studiesSub$study_name
summary(studiesSub$number_of_deployed_locations)

# Filter out ICARUS studies (very few ACC observations)
studiesSub <- studiesSub[grep("ICARUS", studiesSub$study_name, invert=T, value=),]

# Remove duplicated studies and tests
# StudiesSub contains all info per species (therefore the same study could be duplicated if it includes multiple species)
StudiesSub_noDup <- studiesSub[!duplicated(studiesSub[,c("study_name","study_id","timestamp_first_deployed_location","timestamp_last_deployed_location")]),]
# Make sure each study has a unique ID
nrow(StudiesSub_noDup)==length(unique(StudiesSub_noDup$study_id))

# Some studies might be only for testing, we remove them based on the name
if(length(grep("test|Test", StudiesSub_noDup$study_name, value=T))>0){
  StudiesSub_noDup <- StudiesSub_noDup[-grep("test|Test", StudiesSub_noDup$study_name),]
}

write.csv(StudiesSub_noDup, "DataAvailable/GpsAcc_MovebankDatasets/newStudies_toProcess.csv", row.names=F)

#__________________________________________
# Download GPS data from these studies ####

library(plyr)
library(tools)
library(doParallel)
library(data.table)
library(move)
detectCores()
doParallel::registerDoParallel(5)

dir.create("DataDownloaded")

creds <- movebankLogin()

# Import the body mass infos and check that all species names match:
speciesMassInfos <- read.csv("DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy_gps+radar.csv", as.is=T)

# Define function to find errors in the list after the "try"
is.error <- function(x) inherits(x, "try-error")

# Create temporary folder to store the csv data downloaded using the system call
tmpfld <- tempdir()

# Download data per individual so that we can already filter out individuals that don't have acc information
results <- llply(1:nrow(StudiesSub_noDup), function(i) try({
  stRow <- StudiesSub_noDup[i,]
  print(paste0(stRow$study_id," - ",stRow$study_name))
  studyId <- as.numeric(stRow$study_id)
  ## getting license terms of study
  system(paste0('curl -v -u ', paste0(as.vector(creds$headers),collapse=":"), ' -c ./cookies.txt -o ./license_terms.txt "https://www.movebank.org/movebank/service/direct-read?entity_type=event&study_id=', studyId, '"'))
  ## download data using licence. After running this line once, we can use the regular move functions for download. --- CSV files name is the studyID. details about attributes, sensors, individuals, etc can/should be added
  system(paste0('curl -v -u ', paste0(as.vector(creds$headers), collapse=":"), ' -b ./cookies.txt -o ',tmpfld,'/',paste0(studyId,".csv"),' "https://www.movebank.org/movebank/service/direct-read?entity_type=event&study_id=', studyId, '&license-md5=', md5sum("./license_terms.txt"), '"'))
  ## Now the licence has been "accepted" and we can normally downloaad the rest of the information using the move functions
  allInds <- getMovebank("individual", login=creds, study_id=studyId)
  # Exclude individuals that have no ACC data
  allInds_acc <- allInds[grep("GPS.*Acceleration|GPS.*Accessory|acceleration", allInds$sensor_type_ids),]
  # Exclude potential testing individuals
  allInds_acc <- allInds_acc[which(!allInds_acc$local_identifier %in% grep("test|Test", allInds_acc$local_identifier, value=T)),]
  # For the remaining individuals, download the data one by one
  if(nrow(allInds_acc) > 0){
    indNames <- as.character(allInds_acc$local_identifier[!allInds_acc$local_identifier %in% c(NA,"")])
    gps_ls <- lapply(indNames, function(ind)try({
      print(ind)
      # Download gps and check if data are from ornitela or eobs tags
      gps <- getMovebankLocationData(study=studyId, animalName=ind, 
                                     sensorID="GPS", login=creds)
      if(nrow(gps)>0){
        if(length(grep("orn|ornitela", names(gps)))>0){
          gps$deviceType <- "ornitela"
          # In ornitela tags sometimes data are duplicated because sent both via SMS and GPRS. Only GPRS have complete infos so we exclude SMS and keep GPRS (or missing)
          gps <- gps[gps$ornitela.transmission.protocol %in% c(NA,"GPRS"),]
        }else if(length(grep("eobs", names(gps)))>0){ # if eobs filter by status
          gps$deviceType <- "eobs"
          # In eobs tags duplicated are often due to a status different from active, so we only keep observations = A (active) or missing
          if(length(gps$eobs.status)>0){
            gps <- gps[gps$eobs.status %in% c(NA,"","A"),]
            dups <- duplicated(gps[,c("timestamp", "location.long", "location.lat")])
            dupRows <- gps[gps$timestamp %in% gps[dups,"timestamp"],]
            if(nrow(dupRows)>0){
              eventsToDrop <- dupRows$event.id[!dupRows$eobs.status=="A"]
              gps <- gps[!gps$event.id %in% eventsToDrop,]
            }}
        }else if(length(grep("eobs", names(gps)))==0){ #& length(grep("acceleration", names(gps)))==0){ # if not eobs not sure about the filtering
          gps$deviceType <- "other"
        }
        # in all three cases (ornitela, eobs or other tag) add indiv info and body mass
        if(nrow(gps)>0){
          gps$individual.tag.id <- paste0(gps$individual.id,"_",gps$tag.id)
          # Some individuals have no species specified, set that to NA rather than ""
          gps$individual.taxon.canonical.name[gps$individual.taxon.canonical.name==""] <- NA
          # Different individuals in the study might belong to different species, if existing and matching, take the corresponding body mass value from the original study list
          # otherwise body mass value is set to NA
          if(length(unique(gps$individual.taxon.canonical.name))==1 & 
             unique(gps$individual.taxon.canonical.name) %in% speciesMassInfos$matchingSpeciesName){
            gps$BodyMass_value <- unique(speciesMassInfos$BodyMass_value[which(speciesMassInfos$matchingSpeciesName == unique(gps$individual.taxon.canonical.name))])
          }else{gps$BodyMass_value <- NA}
          # Remove NAs from timestamp and coords
          gps <- gps[complete.cases(gps[,c("timestamp","location.long","location.lat")]),]
          if(nrow(gps)>0){return(gps)}}
      }
    }))
    # Remove potential individuals that returned errors during download
    gps_ls <- gps_ls[!vapply(gps_ls, is.error, logical(1))]
    # Exclude empty elements from the list
    gps_ls <- gps_ls[which(!sapply(gps_ls, is.null))]
    # Split the list in two sublists depending on tag type and rbind them and save them separately
    #(some studies have mixed tag types therefore cannot be binded together because of different column names)
    if(length(gps_ls)>0){
      gps_ls_tagType <- split(gps_ls, sapply(lapply(gps_ls, "[",, "deviceType"), unique))
      lapply(names(gps_ls_tagType), function(tag){
        gpsDf <- rbindlist(gps_ls_tagType[[tag]])
        if(tag=="ornitela"){save(gpsDf, file=paste0("DataDownloaded/studyId_",studyId,"_ornitela_movebankDownload_onlyGps.rdata"))
        }else if(tag=="eobs"){save(gpsDf, file=paste0("DataDownloaded/studyId_",studyId,"_eobs_movebankDownload_onlyGps.rdata"))
        }else if(tag=="other"){save(gpsDf, file=paste0("DataDownloaded/studyId_",studyId,"_other_movebankDownload_onlyGps.rdata"))}
      })
    }
  }
}), .parallel=F) #working in parallel doesn't work for the download

# Check errors
# is.error <- function(x) inherits(x, "try-error")
# succeeded <- !vapply(results, is.error, logical(1))
# StudiesSubSub_toDo$study_name[which(StudiesSubSub_toDo!=succeeded)]

#______________________________________________________________________
# Explore/find common sampling frequency to subset all datasets ####

# Make some timelag plots to decide on the sampling interval we want to use before associating the acc values
dir.create("Plots")

fls <- list.files("DataDownloaded", "eobs.*onlyGps|other.*onlyGps", full.names=T) # excluding ornitela

pdf("Plots/timeLag_allStudies_gps.pdf")
TLsumm_ls <- lapply(fls, function(f){
  print(f)
  load(f) # object gpdDf
  TL <- unlist(sapply(split(gpsDf, gpsDf$individual.tag.id), function(ind){
    return(as.numeric(difftime(ind$timestamp[-1],ind$timestamp[-nrow(ind)], units="mins")))
  }))
  hist(TL, breaks="FD", main=paste0(unique(gpsDf$study.name),"\n median TL = ",median(TL)," min"), xlim=c(0,2*60), xlab="Timelag in minutes (TL)")
  return(summary(TL))
})
dev.off()


#____________________________________________
# Download ACC data for same studies ####

options(digits.secs=4)

library(move)
library(data.table)
library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(3)
is.error <- function(x) inherits(x, "try-error")

library(RCurl)
curlSetOpt(timeout=400) #to solve timeout problem

fls <- list.files("DataDownloaded", "eobs.*onlyGps|other.*onlyGps", full.names=T) # excluding ornitela because acc data have to be downloaded directly from ornitela website

creds <- movebankLogin()

results <- llply(fls, function(f)try({
  # Extract study name from file names of downloaded gps files
  studyId <- as.numeric(strsplit(f, "_")[[1]][2]) 
  deviceType <- strsplit(f, "_")[[1]][3]
  print(f)
  # List names of the individuals included in the study and filter only those having ACC data
  allInds <- getMovebank("individual", login=creds, study_id=studyId)
  allInds_acc <- allInds[grep("GPS.*Acceleration|GPS.*Accessory|acce.*gps", allInds$sensor_type_ids),]
  allInds_acc <- allInds_acc[which(!allInds_acc$local_identifier %in% grep("test|Test", allInds_acc$local_identifier, value=T)),]
  # For the remaining individuals, download the data one by one
  if(nrow(allInds_acc) > 0){
    indInfo <- data.frame(indName=as.character(allInds_acc$local_identifier[!allInds_acc$local_identifier%in%c(NA,"")]),
                          indId=as.numeric(allInds_acc$id[!allInds_acc$local_identifier%in%c(NA,"")]))
    accLs <- lapply(indInfo$indId, function(ind)try({  #sample(indInfo$indId, 25) - if file is too big sample few individuals
      print(ind)
      # Download only the timestamps of the acc data
      Trange <- getMovebank(entity="event", study_id=studyId, individual_id=ind, sensor_type_id=2365683, login=creds,
                            attributes=c("timestamp"))#, timestamp_end="20220610235959") #up to 10th of June 2022
      Trange <- as.vector(Trange$timestamp)
      if(length(Trange)>0){
        # Split the timestamps in chunks of 100k lines per download (otherwise the download breaks)
        dlsq <- cbind(seq(0, length(Trange), by=100000)+1, c(seq(0, length(Trange), by=100000), length(Trange))[-1])
        # Create start and end timestamps for download
        DLDtimeSeq <- data.frame(startTime=Trange[dlsq[,1]], endTime=Trange[dlsq[,2]])
        # Download the chunks one by one
        ACCdwld <- rbindlist(lapply(1:nrow(DLDtimeSeq), function(j){
          chunkACC <- getMovebankNonLocationData(study=studyId, animalName=ind, sensorID="Acceleration", login=creds, 
                                                 timestamp_start=gsub(".", "", gsub("-|:| ", "", DLDtimeSeq$startTime[j]), fixed=T), 
                                                 timestamp_end=gsub(".", "", gsub("-|:| ", "", DLDtimeSeq$endTime[j]), fixed=T))
          chunkACC$timestamp <- as.POSIXct(strptime(chunkACC$timestamp, "%Y-%m-%d %H:%M:%OS", tz="UTC"))
          return(chunkACC)
        }))
        return(as.data.frame(ACCdwld))
      }
    }))
     # Remove potential individuals that returned errors during download
    accLs <- accLs[!vapply(accLs, is.error, logical(1))]
    # Exclude empty elements from the list, bind and save
    accLs <- accLs[which(!sapply(accLs, is.null))]
    if(length(accLs)>0){
      accDf <- rbindlist(accLs)
      save(accDf, file=paste0("DataDownloaded/studyId_",studyId,"_",deviceType,"_movebankDownload_onlyAcc.rdata"))
    }
  }
}), .parallel=F) #parallel gives errors with download


#____________________________________________________
# Associate GPS and ACC and calculate mean vedba ####

#____________
# Download the tag reference table for all studies to make sure to know how to treat the different ACC data formats
library(move)
library(data.table)

setwd("...")

fls_gps <- list.files("DataDownloaded", "onlyGps", full.names=T)
fls_acc <- list.files("DataDownloaded", "onlyAcc", full.names=T)
# Make sure to select only GPS files (of each individuals) for which ACC files also exist
fls <- fls_gps[sapply(strsplit(fls_gps, "_"),"[",2) %in% sapply(strsplit(fls_acc, "_"),"[",2)]

accStudyIds <- as.numeric(unique(sapply(strsplit(fls, "_"), "[", 2)))
tagRef_df <- as.data.frame(rbindlist(lapply(accStudyIds, function(studyId){
  tagDf <- getMovebank("tag", login=creds, study_id=studyId)
  tagDf <- tagDf[,c("local_identifier","manufacturer_name","weight")]
  tagDf$studyId <- studyId
  names(tagDf)[1] <- "tag_local_identifier"
  return(tagDf)
})))

tagRef_df$tag_local_identifier <- as.character(tagRef_df$tag_local_identifier)
tagRef_df$studyId <- as.character(tagRef_df$studyId)
tagRef_df$manufacturer_name <- as.character(tagRef_df$manufacturer_name)
tagRef_df$manufacturer_name[tagRef_df$manufacturer_name==""] <- NA
tagRef_df$manufacturer_name[grep(tagRef_df$manufacturer_name, pattern="obs|Obs|OBS|E-ops|eObd GmBH")] <- "eobs"
tagRef_df$manufacturer_name[grep(tagRef_df$manufacturer_name, pattern="UKn|Uni KN|Konstanz|Uni Kn")] <- "UniKonstanz"
tagRef_df$manufacturer_name[grep(tagRef_df$manufacturer_name, pattern="Microwave")] <- "Microwave"
table(tagRef_df$manufacturer_name)
tagRef_ls <- split(tagRef_df, tagRef_df$studyId)
names(tagRef_ls)
sapply(tagRef_ls, function(t)unique(as.character(t$manufacturer_name)))
lapply(tagRef_ls, head)
save(tagRef_ls, file="DataDownloaded/referenceTable_manufacturer_allStudies.rdata")

#____________
#Now associate gps and acc variables (obtained depending on the tag manufacturer)

library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(5)
library(move)
library(lubridate)
library(data.table)

source("Scripts/COT_publ/COT_functions.R") # load own function library (ACCtoGPS function)
is.error <- function(x) inherits(x, "try-error")

options(digits=6) #for decimals


load("DataDownloaded/referenceTable_manufacturer_allStudies.rdata") #obj tagRef_ls
fls_gps <- list.files("DataDownloaded", "onlyGps", full.names=T)
fls_acc <- list.files("DataDownloaded", "onlyAcc", full.names=T)
# Select only GPS files for which ACC files also exist
fls <- fls_gps[sapply(strsplit(fls_gps, "_"),"[",2) %in% sapply(strsplit(fls_acc, "_"),"[",2)]

# check the remaining files yet to be processed
flsDone <- list.files("DataDownloaded", "gpsAcc.rdata", full.names=T)
if(length(flsDone)>0){
  flsToDo <- fls[!fls %in% grep(paste(sapply(strsplit(flsDone, "_"),"[",2),collapse="|"), fls, value=T)] # subset of fls that still haven't gone through this step
}else{flsToDo <- fls}

results <- lapply(flsToDo, function(f)try({
  studyId <- strsplit(f, "_")[[1]][2]
  # Import accelerometry data downloaded in the previous step
  load(grep(studyId, fls_acc, value=T)) #data.frame, object accDf
  # Use the tag reference df to assign device type to each individual
  studyTagRef <- tagRef_ls[[studyId]][tagRef_ls[[studyId]]$tag_local_identifier %in% unique(accDf$tag_local_identifier),]
  deviceType <- studyTagRef$manufacturer_name[!is.na(studyTagRef$manufacturer_name)]
  if(length(deviceType)==0){deviceType <- strsplit(f, "_")[[1]][3]}
  deviceType <- paste(unique(deviceType), collapse="|")
  # Import the gps data for that study and split by individual
  load(f) #data.frame, object gpsDf
  print(paste0(unique(gpsDf$study.name)," - ",studyId))
  # Update deviceType with info present in the tag reference table by merging the tag ref table with the gps df so that each tag is allowed to have a different device type
  gpsDf$tag.local.identifier <- as.character(gpsDf$tag.local.identifier)
  gpsDf <- merge(gpsDf, studyTagRef[,c("tag_local_identifier","manufacturer_name")], by.x="tag.local.identifier", by.y="tag_local_identifier", all.x=T)
  gpsDf$deviceType <- as.character(gpsDf$manufacturer_name)
  if(all(is.na(gpsDf$deviceType))){gpsDf$deviceType <- deviceType}
  if(length(unique(studyTagRef$manufacturer_name[!is.na(studyTagRef$manufacturer_name)]))==1){gpsDf$deviceType[which(is.na(gpsDf$deviceType))] <- deviceType}
  gpsDf$manufacturer_name <- NULL
  #Work on each individual separately
  gps_ls <- split(gpsDf, as.character(gpsDf$individual.local.identifier))
  acc_ls <- split(accDf, as.character(accDf$individual_local_identifier))
  # Remove the DFs (now in lists) to free up memory
  rm(gpsDf, accDf)
  # Work on each individual, but only those that are in both the GPS and the ACC lists
  gps_ls <- gps_ls[names(gps_ls) %in% names(acc_ls)]
  gpsAcc_ls <- lapply(names(gps_ls), function(ind)try({
    print(ind)
    # Extract acc and gps for each individual and transform data.table in data.frame
    gps <- as.data.frame(gps_ls[[ind]])
    acc <- as.data.frame(acc_ls[[ind]])
    # Subsample gps to the minimum sampling frequency we decided (5 min)
    gps <- gps[!duplicated(round_date(gps$timestamp,"5 mins")),]
    # Order both datasets by timestamp
    acc <- acc[order(acc$timestamp),]
    gps <- gps[order(gps$timestamp),]
    # Subset only acc locations within the time range of the gps for the interpolation
    acc <- acc[acc$timestamp > min(gps$timestamp) & acc$timestamp < max(gps$timestamp),]
    # Change event_id col in acc_event_id not to confuse them between gps and acc
    names(acc)[names(acc) == "event_id"] <- "acc_event_id"
    acc$acc_event_id <- as.character(acc$acc_event_id)
    # nrow(acc)==length(unique(acc$acc_event_id))
    if(nrow(acc) > 0){
      # This procedure is for "eobs-like" data (all acc value together in one character string)
      if(deviceType %in% c("eobs","madebytheo","milsar") & unique(gps$deviceType) %in% c("eobs","madebytheo","milsar","")){ 
      # Extract column names (sometimes they differ depending on tags)
      if(length(grep("accelerations_raw", names(acc)))>1){
        accRawCol <- grep("eobs_accelerations_raw", names(acc), value=T)
        axesCol <- grep("eobs_acceleration_axes", names(acc), value=T)
      }else{
        accRawCol <- grep("accelerations_raw", names(acc), value=T)
        axesCol <- grep("acceleration_axes", names(acc), value=T)
      }
      if(median(as.numeric(nchar(as.character(acc[,accRawCol]))))>100){
        # If there are duplicated timestamps, remove those that are not "Active" (they contain less information)
        dups <- duplicated(gps[,c("timestamp", "location.long", "location.lat")])
        dupRows <- gps[gps$timestamp %in% gps[dups,"timestamp"],]
        eventsToDrop <- dupRows$event.id[!dupRows$eobs.status=="A"]
        gps <- gps[!gps$event.id %in% eventsToDrop,]
        # Exclude ACC data that have only X or Y axis
        acc <- acc[which(!as.character(acc[, axesCol]) %in% c("X","Y")),]
        # Exclude rows with missing acc data
        acc <- acc[which(complete.cases(acc[, accRawCol])),]
        # Associate to each gps location the ACC EVENT ID of the acc values closest in time (within 5 min) using self made function
        timeTolerance <- 5*60 # 5 mins in seconds
        gpsAcc <- ACCtoGPS(acc, gps, timeTolerance) #table(is.na(gpsAcc$acc_event_id))
        # Subset the acc dataset only to the acc_event_id associated to the gps information
        accSub <- acc[which(acc$acc_event_id %in% gpsAcc$acc_event_id),]
        # Calculate acc infos (vedba etc) for this subset of acc data (split acc strings and calculate mean and cumulative vedba)
        if(nrow(accSub)>0){ 
          #if(length(unique(accSub[, axesCol]))>1 ){ #Continue only if number of acc axes doesn't vary within the same individual
          accDf_vedba <- createVedbaDF_eobs(acc=accSub, accEventCol="acc_event_id", accRawCol=accRawCol, axesCol=axesCol)
          # Merge the average vedba values to the gps based on acc_event_id
          gpsAcc <- merge(gpsAcc, accDf_vedba[,c("acc_event_id","n_samples_per_axis","acc_burst_duration_s","acc_sampl_freq","meanVedba","cumVedba")], 
                          by="acc_event_id", all.x=T)
          #nrow(gpsAcc[which(gpsAcc$meanVedba==0),])
          #print(table(is.na(gpsAcc$meanVedba)))
        }
      }
    }
      if(nrow(gpsAcc[which(gpsAcc$meanVedba==0),])>0){warning("This study contains meanVedba values = 0. Please double check!")}
      return(gpsAcc)
      }
  }))
  # Remove potential individuals that returned errors during download
  gpsAcc_ls <- gpsAcc_ls[!vapply(gpsAcc_ls, is.error, logical(1))]
  # Exclude empty elements from the list, bind and save
  gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
  if(length(gpsAcc_ls)>0){
    gpsDf <- rbindlist(gpsAcc_ls)
    save(gpsDf, file=paste0("DataDownloaded/studyId_",studyId,"_",deviceType,"_movebankDownload_gpsAcc.rdata"))
  }
}))


#__________________________________________________________________________________
# ONLY FOR EOBS: re-calculate mean vedba using ACC data tranformed in Gs ####

# We need to repeat this step for Eobs devices in order to re-calculate VeDBA using the transformed ACC values (in Gs) rather than the raw values.
library(moveACC)
# We will use the gpsAcc eobs files that were already associated in the previous step, but we will add two additional columns with the transformed meanVedba and transformed cumVedba
fls_acc <- list.files("DataDownloaded", "onlyAcc", full.names=T)
fls_gps <- list.files("DataDownloaded", "gpsAcc.rdata", full.names=T)
flsToTransform <- grep("eobs", fls_gps, value=T)
flsDone <- list.files("DataDownloaded", "gpsAcc_transf", full.names=T)
flsToTransform <- flsToTransform[!sapply(strsplit(flsToTransform, "_"),"[",2) %in% sapply(strsplit(flsDone, "_"),"[",2)]
#flsToTransform <- flsToTransform[sapply(strsplit(flsToTransform, "_"),"[",2) %in% sapply(strsplit(fls_acc, "_"),"[",2)]

lapply(flsToTransform, function(f)try({
  print(f)
  studyId <- strsplit(f, "_")[[1]][2]
  # Import gps and accelerometry data downloaded in the previous step
  load(f) #object gpsDf
  load(grep(studyId, fls_acc, value=T)) #data.frame, object accDf
  # transform data.table in data.frame
  gpsDf <- as.data.frame(gpsDf)
  accDf <- as.data.frame(accDf)
  # extract column with individual id
  if("individual_local_identifier" %in% names(accDf)){indCol <- "individual_local_identifier"}else{indCol <- "local_identifier"}
  #Work on each individual separately
  gps_ls <- split(gpsDf, as.character(gpsDf$individual.local.identifier))
  acc_ls <- split(accDf, as.character(accDf[,indCol])) #because it's a data.table object
  # Remove the DFs (now in lists) to free up memory
  rm(gpsDf, accDf)
  # Work on each individual, but only those that are in both the GPS and the ACC lists
  gps_ls <- gps_ls[names(gps_ls) %in% names(acc_ls)]
  acc_ls <- acc_ls[names(acc_ls) %in% names(gps_ls)]
  gpsAcc_ls <- llply(names(gps_ls), function(ind)try({
    print(ind)
    # Extract acc and gps for each individual and transform data.table in data.frame
    gps <- gps_ls[[ind]]
    acc <- acc_ls[[ind]]
    # Order both datasets by timestamp
    acc <- acc[order(acc$timestamp),]
    gps <- gps[order(gps$timestamp),]
    # Subset only acc locations within the time range of the gps for the interpolation
    acc <- acc[acc$timestamp > min(gps$timestamp) & acc$timestamp < max(gps$timestamp),]
    # Change event_id col in acc_event_id not to confuse them between gps and acc
    names(acc)[names(acc) == "event_id"] <- "acc_event_id"
    acc$acc_event_id <- as.character(acc$acc_event_id)
    if(nrow(acc) > 0){
      # extract the column name of the raw acc values and axes and tag id
        accRawCol <- grep("accelerations_raw|accelerations.raw", names(acc), value=T)
        axesCol <- grep("acceleration.axes|acceleration_axes", names(acc), value=T)
        tagidCol <- grep("tag.local.identifier|tag_local_identifier", names(acc), value=T)
        # make sure the tag id column is numeric (but careful with factors)
        acc[,tagidCol] <- as.integer(as.character(acc[,tagidCol]))
        # if more than one column exists for acc raw and axes, remove the ones with more NAs (otherwise the function will detect multiple columns with raw acc)
        if(length(accRawCol)==2){
          rawToDrop <- accRawCol[which.max(c(length(which(is.na(acc[,accRawCol[1]]))), length(which(is.na(acc[,accRawCol[2]])))))] #compare the number of NAs and exclude the one with more
          acc[,rawToDrop] <- NULL
          accRawCol <- grep("accelerations_raw|accelerations.raw", names(acc), value=T)
        }
        if(length(axesCol)==2){
          axesToDrop <- axesCol[which.max(c(length(which(is.na(acc[,axesCol[1]]))), length(which(is.na(acc[,axesCol[2]])))))] #compare the number of NAs and exclude the one with more
          acc[,axesToDrop] <- NULL
          axesCol <- grep("acceleration.axes|acceleration_axes", names(acc), value=T)
        }
        if(median(as.numeric(nchar(as.character(acc[, accRawCol]))))>100){
        # Exclude ACC data that have only X or Y axis
        acc <- acc[which(!as.character(acc[, axesCol]) %in% c("X","Y")),]
        # Exclude rows with missing acc data
        acc <- acc[which(complete.cases(acc[, accRawCol])),]
        # !!ACC were already associated to the GPS in the previous step, so we can simply use the same acc_event_id to retrieve the correct ACC burst 
        # Subset the acc dataset only to the acc_event_id associated to the gps information
        accSub <- acc[which(acc$acc_event_id %in% gps$acc_event_id),]
        # Transform raw eobs values in Gs (it assumes the column containing raw values is names as in accRawCol defined above)
        # for eobs tags with numbers smaller than 2241 a sensitivity setting 'low' or 'high' has to be provided.
        if(any(accSub[,tagidCol] < 2241)){sens <- data.frame(TagID=unique(accSub[accSub[,tagidCol] < 2241, tagidCol]), sensitivity="low")}else{sens <- NULL}
        accTransf <- TransformRawACC(accSub, sensitivity.settings = sens, units="g", showProgress=F)
        # Calculate acc infos (vedba etc) for this subset of acc data (split acc strings and calculate mean and cumulative vedba)
        if(nrow(accTransf)>0){
          # we recalculate vedba, in this case using the transformed acc values, so vedba will be in Gs
          accDf_vedba <- createVedbaDF_eobs(acc=accTransf, accEventCol="acc_event_id", accRawCol="accelerationTransformed", axesCol=axesCol)
          names(accDf_vedba)[names(accDf_vedba) %in% c("meanVedba","cumVedba")] <- c("meanVedba_Gs","cumVedba_Gs")
          # Merge the average vedba values to the gps based on acc_event_id
          gpsAcc <- merge(gps, accDf_vedba[,c("acc_event_id","meanVedba_Gs","cumVedba_Gs")], 
                          by="acc_event_id", all.x=T)
          #nrow(gpsAcc[which(gpsAcc$meanVedba==0),])
          #print(table(is.na(gpsAcc$meanVedba)))
        }
      }
      if(nrow(gpsAcc[which(gpsAcc$meanVedba==0),])>0){warning("This study contains meanVedba values = 0. Please double check!")}
      return(gpsAcc)
    }
    }), .parallel=F)
  # Remove potential individuals that returned errors during download
  gpsAcc_ls <- gpsAcc_ls[!vapply(gpsAcc_ls, is.error, logical(1))]
  # Exclude empty elements from the list, bind and save
  gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
  if(length(gpsAcc_ls)>0){
    gpsDf <- rbindlist(gpsAcc_ls)
    save(gpsDf, file=paste0("DataDownloaded/studyId_",studyId,"_eobs_movebankDownload_gpsAcc_transf.rdata"))
  }
}))


