
#_________________________________________________________
# Transform in move, calculate track geometry metrics ####
#_________________________________________________________


library(move)
library(sf)
library(amt)
library(data.table)
library(plyr)
library(doParallel)

library(lubridate)
options(digits=6)

setwd("...")
source("Scripts/COT_publ/COT_functions.R") # load own function library

flsEobs <- list.files("DataDownloaded", "gpsAcc_transf", full.names=T)
flsOthers <- grep("eobs", list.files("DataDownloaded", "gpsAcc", full.names=T), invert = T, value=T)
fls <- c(flsEobs, flsOthers)

dir.create("DataProcessed")

flsDone <- list.files("DataProcessed", pattern="dfGpsAcc_geom_10min", full.names = T) # files that already went through this step
if(length(flsDone)>0){
  flsToDo <- fls[!fls %in% grep(paste(sapply(strsplit(flsDone, "_"),"[",2),collapse="|"), fls, value=T)] # subset of fls that still haven't gone through this step
}else{flsToDo <- fls}

doParallel::registerDoParallel(2)
results <- llply(flsToDo, function(f)try({
  load(f) #data.frame object gpsDf
  if(any(class(gpsDf)=="data.table")){gpsDf <- as.data.frame(gpsDf)}  #Transform data.table into data.frame objects
  studyId <- as.numeric(strsplit(f, "_")[[1]][2])
  print(paste0(unique(gpsDf$study.name)," - ",studyId))
  gps_ls <- split(gpsDf, as.character(gpsDf$individual.tag.id))
  gpsAcc_geom_ls <- lapply(gps_ls, function(ind){
    print(unique(ind$individual.tag.id))
    # Thin the data to 10+-5 min
    # (for eobs this step was already done in previous step before associating acc)
    ind <- ind[!duplicated(ind[,c("timestamp","location.long","location.lat","meanVedba_Gs")]),] #in some cases some duplicates are left, let's remove them
    indAmt <- mk_track(tbl=ind, all_cols=T,
                       .x=location.long, .y=location.lat, crs = st_crs(4326),
                       .t=timestamp, order_by_ts = T, check_duplicates = T)
    indAmt <- track_resample(indAmt, rate = minutes(10), tolerance = minutes(5), start = 1)
    m <- as_move(indAmt)
    # Remove outliers based on speed (max 50 m/s)
    while(any(move::speed(m) >= 50) == T){
      m <- m[which(c(NA, move::speed(m)) < 50)]
    }
    # # Check that there is no duplicates before transforming into a move object
    # dup_ls <- getDuplicatedTimestamps(x=as.factor(ind$individual.tag.id),
    #                                   timestamps=ind$timestamp)
    # if(length(dup_ls)==0){ 
    #   # Transform in move object
    #   m <- move(x=ind$location.long, y=ind$location.lat, proj=CRS("+proj=longlat +ellps=WGS84"),
    #             time=ind$timestamp, data=ind)
    if(n.locs(m)>0){
      # Extract altitude column
      heightCol <- grep("height_above|height.above", names(m@data), value=T)
      height <- m@data[, heightCol]
      # If there is no altitude column means that there is no altitude information
      if(length(height)==0){
        height <- rep(NA, nrow(m@data))
        # If there are multiple height columns keep the one with the least number of NAs
      }else if(length(heightCol)>1){
        lessNA <- which.min(sapply(1:ncol(height), function(c){
          length(which(is.na(height[,c])))
        }))
        height <- height[,lessNA]
      }
      # Add variables about the track geometry
      m$altitudeDiff <- c(NA, (height[-1] - height[-length(height)]))
      m$timeLag_min <- c(NA, timeLag(m, units="mins"))
      m$vertSpeed_ms <- m$altitudeDiff/(m$timeLag_min*60)
      m$stepLength_m <- c(NA, move::distance(m))
      m$groundSpeed_ms <- c(NA, move::speed(m))
      m$segmentDir <- c(NA, direction360(bearing(m)[-nrow(m)]))
      m$turnAngle <- c(NA, turnAngleGc(m), NA)
      # achieving the same without package move
      #m$stepLength_m <- c(NA, raster::pointDistance(m[-nrow(m),c("location.long","location.lat")], m[-1,c("location.long","location.lat")], longlat=T))
      #m$groundSpeed_ms <- m$stepLength_m/(m$timeLag_min*60)
      #m$turnAngle <- c(m$segmentDir180[-1]-m$segmentDir180[-nrow(m)], NA) #without move you need segmentDir -180/180
      #m$segmentDir <- c(NA, direction360(bearing(m[-nrow(m),c("location.long","location.lat")], m[-1,c("location.long","location.lat")]))) #without move package
      return(as.data.frame(m))
    }#}
  })
  # Exclude empty elements from the list, and bind the dataframes of all the non null individuals
  if(is.null(gpsAcc_geom_ls)==F){
    gpsAccDF_geom <- as.data.frame(rbindlist(gpsAcc_geom_ls[which(!sapply(gpsAcc_geom_ls, is.null))], use.names = T))
    # Keep only data points for which there are ACC information available
    accInfo <- gpsAccDF_geom[, grep("meanVedba|cumAccAxes", names(gpsAccDF_geom), value=T)]
    gpsAccDF_geom <- gpsAccDF_geom[complete.cases(accInfo),]
    if(nrow(gpsAccDF_geom)>0){
      # Split and save per species in case there is more than one species in the same study
      spCol <- grep("taxon_canonical_name|taxon.canonical.name", names(gpsAccDF_geom), value=T)
      lapply(split(gpsAccDF_geom, as.character(gpsAccDF_geom[,spCol])), function(df_geom){
        spName <- paste(strsplit(unique(df_geom[,spCol])," ")[[1]],collapse="-")
        save(df_geom, file=paste0("DataProcessed/studyId_",studyId,"_",spName,"_",unique(df_geom$deviceType),"_dfGpsAcc_geom_10min.RData"))
      })
    }
  }
}), .parallel=T)

lapply(results, head, 1) #check errors
