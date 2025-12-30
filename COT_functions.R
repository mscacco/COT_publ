
#___________________________________________
## ALL FUNCTIONS NEEDED FOR THE COT PROJECT

#_____________________________
#Function to load a rdata file with the name you want (so that you don't have to load it/save it with another name/cancel the one you imported)
load.rdata <- function(file.path){
  #loads an RData file, and returns it
  load(file.path)
  get(ls()[ls() != "file.path"]) #return everything except the object called with the name it was saved with 
}

#___________________________________________
## Function to assign each acc point to the gps info nearest in time - Modified from Anne's version (associate_ACC&GPS_Anne) #####
#E.g. to add lat-long and altitude of the gps to the closest acc

#function arguments are acc data, gps data, columns of gps that you want to associate to the acc,
#ColsToCreate = name of the new columns to create in the acc dataset (optional, by default are the same names of the gps ColsToAssociate)
#and a time tolerance to associate them in seconds (e.g. within 30 secs before or after my acc point)

## Needs library(plyr) and library(doParallel) to run because of the llply

#function arguments are acc data, gps data, columns of acc that you want to associate to the gps
#and a time tolerance to associate them in seconds (e.g. within 30 secs before or after my acc point)
#Here GPS is the reference, the function returns the GPS data with associated, when available, the ID of the closest ACC event
ACCtoGPS <- function(ACCdata,GPSdata,
                     timeTolerance,
                     accEventCol=NULL,
                     ColsToAssociate=NULL,
                     ACCtimeCol="timestamp", GPStimeCol="timestamp"){
  require(data.table)
  if(is.null(accEventCol)){accEventCol <- grep("event_id|event.id", names(ACCdata), value=T)}
  if(length(accEventCol)==0){
    ACCdata$event_id <- 1:nrow(ACCdata)
    accEventCol <- "event_id"}
  if(all(class(ACCdata[,ACCtimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F) | 
     all(class(GPSdata[,GPStimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F)){
    stop("ACC or GPS time column is not in POSIXct format.")}
  #create the empty columns that I want to fill in during the loop
  GPSdata$acc_event_id <- NA
  GPSdata$diff_acc_time_s <- NA
  GPSdata$acc_closest_timestamp <- NA
  ACCdata <- ACCdata[ACCdata[,ACCtimeCol] > min(GPSdata[,GPStimeCol]) & ACCdata[,ACCtimeCol] < max(GPSdata[,GPStimeCol]),]
  if(nrow(ACCdata)==0){stop("There are no ACC data available in the GPS time range. Consider chacking that the two datasets are in the same time zone.")}
  GPSdata <- as.data.frame(rbindlist(lapply(1:nrow(GPSdata), function(h){
    #create a subset of acc that occurr within a certain time interval from each gps point (+/- timeTolerance)
    gps.time <- GPSdata[h,GPStimeCol]
    acc.sub <- ACCdata[ACCdata[,ACCtimeCol] > gps.time-timeTolerance & ACCdata[,ACCtimeCol] < gps.time+timeTolerance,]
    #there could be no acc data in that close interval, but if there are some (nrow > 1) then we can associate them to the gps info:
    if(nrow(acc.sub) >= 1){
      timeDiff <- abs(difftime(acc.sub[,ACCtimeCol], GPSdata[h,GPStimeCol], units="secs")) #calculates the time difference
      min.diff <- which.min(timeDiff) #selects the row of the nearest point in time of the acc to each gps point (h)
      ## take the acc event id from the point corresponding to the minimum time difference and associate it to the gps point h
      GPSdata$acc_event_id[h] <-  acc.sub[min.diff, accEventCol]
      GPSdata$diff_acc_time_s[h] <- round(min(as.difftime(timeDiff, units="secs")), digits=3)
      GPSdata$acc_closest_timestamp[h] <- as.character(acc.sub[min.diff,ACCtimeCol])
    }
    return(GPSdata[h,])
  })))
  if(all(is.na(GPSdata$acc_closest_timestamp))){warning("No ACC data to associate to any of the GPS burst given this time tolerance. All associated ACC information set to NA. Consider increasing the time tolerance.")}
  GPSdata$acc_closest_timestamp <- as.POSIXct(GPSdata$acc_closest_timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  if(!is.null(ColsToAssociate)){
    GPSdata <- merge(GPSdata, ACCdata[,c(accEventCol,ColsToAssociate)], by.x="acc_event_id", by.y=accEventCol, all.x=T)
    GPSdata <- GPSdata[order(GPSdata[,GPStimeCol]),]
  }
  return(as.data.frame(GPSdata))  #the function returns the GPS dataset with the additional ACC columns + a time.diff column
}

#_____________________________
# Function to calculate vedba per acc observation, starting from a X Y and Z acc axis
# static acceleration is not based on a rolling mean but is calculated as a mean per burst

calculateVedba <- function(accX, accY, accZ=NULL){
  if(is.null(accX)|is.null(accY)){stop("Less than 2 axes available: VeDBA cannot be computed.")}
  if(is.null(accZ)){
    vedba <- sqrt((accX-mean(accX))^2 + (accY-mean(accY))^2)
    warning("Z axis is missing: VeDBA is being calculated on 2 axes.")
  }else{vedba <- sqrt((accX-mean(accX, na.rm=T))^2 + (accY-mean(accY, na.rm=T))^2 + (accZ-mean(accZ, na.rm=T))^2)}
  return(vedba)
}
# This variation of the same function allows for calculating DBA on the Z axis (the one above excludes data were only the Z axis is available)
calculateVedba_alsoZ <- function(accX=NULL, accY=NULL, accZ=NULL, nAxes=3){
  if(nAxes==1 & is.null(accZ)){
    stop("Only 1 axis available, which is not the Z axis: VeDBA cannot be computed.")
  }else if(nAxes==1 & !is.null(accZ)){
    vedba <- sqrt((accZ-mean(accZ))^2)
    warning("X and Y axes are missing. (Ve)DBA is being calculated only on the Z axis.")
  }else if(nAxes==2 & is.null(accZ)){
    vedba <- sqrt((accX-mean(accX))^2 + (accY-mean(accY))^2)
    warning("Z axis is missing: VeDBA is being calculated on 2 axes.")
  }else{vedba <- sqrt((accX-mean(accX, na.rm=T))^2 + (accY-mean(accY, na.rm=T))^2 + (accZ-mean(accZ, na.rm=T))^2)}
  return(vedba)
}

#_____________________________
# Function to calculate mean vedba per burst for eobs tags (burst data with all axes in one string)

createVedbaDF_eobs <- function(acc, 
                               accEventCol=grep("acc_event_id|acc.event.id", names(acc), value=T), 
                               axesCol=grep("acceleration_axes|acceleration.axes", names(acc), value=T), 
                               accRawCol=grep("accelerations_raw|accelerations.raw", names(acc), value=T), 
                               sampFreqCol=grep("acceleration_sampling_frequency_per_axis|acceleration.sampling.frequency.per.axis", names(acc), value=T)){
  # if multiple columns are available for the same information, take the one with the least amount of NAs
  if(length(axesCol)==2){
    axesCol <- axesCol[which.min(c(length(which(is.na(acc[,axesCol[1]]))), length(which(is.na(acc[,axesCol[2]])))))] #compare the number of NAs and exclude the one with more
  }
  if(length(accRawCol)==2){
    accRawCol <- accRawCol[which.min(c(length(which(is.na(acc[,accRawCol[1]]))), length(which(is.na(acc[,accRawCol[2]])))))] #compare the number of NAs and exclude the one with more
  }
  if(length(sampFreqCol)==2){
    sampFreqCol <- sampFreqCol[which.min(c(length(which(is.na(acc[,sampFreqCol[1]]))), length(which(is.na(acc[,sampFreqCol[2]])))))] #compare the number of NAs and exclude the one with more
  }
  #Continue only if number of acc axes doesn't vary within the same individual
  if(nrow(acc)==0){stop("The acc dataset has 0 observations.")
  }else if(nrow(acc)>0){
    if(length(unique(acc[,axesCol]))>1){
      warning("Note that the ACC observations in this dataset have variable number of axes.")}
    #Create an empty dataframe with one row per burst to fill in with the vedba values
    accDf_vedba <- data.frame(acc_event_id=acc[,accEventCol], 
                              #n_axes=NA,
                              n_samples_per_axis=NA, acc_burst_duration_s=NA, acc_sampl_freq=NA,
                              meanVedba=NA, cumVedba=NA, 
                              stringsAsFactors = F)
    #fill it in with mean and cumulative vedba, number of samples per axis and burst duration
    acc[,axesCol] <- as.character(acc[,axesCol])
    for(j in 1:nrow(acc)){
      Naxes <- nchar(as.character(acc[j, axesCol]))
      accMx <- matrix(as.numeric(unlist(strsplit(as.character(acc[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
      n_samples_per_axis <- nrow(accMx)
      acc_burst_duration_s <- n_samples_per_axis/acc[j, sampFreqCol]
      if(nchar(acc[j, axesCol])==1 & unique(acc[j, axesCol])!="Z"){stop("Only 1 axis available, and it is not Z.")}
      if(nchar(acc[j, axesCol])==1 & unique(acc[j, axesCol])=="Z"){
        vedba <- calculateVedba_alsoZ(accZ=accMx[,1], nAxes=Naxes) #needs to specify that this is axis Z
      }
      if(nchar(acc[j, axesCol])==2){
        vedba <- calculateVedba_alsoZ(accMx[,1], accMx[,2], nAxes=Naxes)
      }
      if(nchar(acc[j, axesCol])==3){
        vedba <- calculateVedba_alsoZ(accMx[,1], accMx[,2], accMx[,3], nAxes=Naxes)
      }
      accDf_vedba[j, c("n_samples_per_axis","acc_burst_duration_s","acc_sampl_freq",
                       "meanVedba","cumVedba")] <- c(n_samples_per_axis, acc_burst_duration_s, acc[j, sampFreqCol],
                                                     mean(vedba, na.rm=T), sum(vedba, nrm=T))
      # Modification including the column number of axes:
      # accDf_vedba[j, c("n_axes","n_samples_per_axis","acc_burst_duration_s","acc_sampl_freq",
      #                  "meanVedba","cumVedba")] <- c(Naxes, n_samples_per_axis, acc_burst_duration_s, acc[j, sampFreqCol],
      #                                                mean(vedba, na.rm=T), sum(vedba, nrm=T))
    }
    return(accDf_vedba)
  }
}

#_____________________________
# Function to calculate mean vedba per burst for uvaBits tags (bursts but axes in separate columns)
# 
# 
# createVedbaDF_uvaBits <- function(acc, 
#                                accEventCol=grep("acc_event_id|acc.event.id", names(acc), value=T), 
#                                axisX=grep("acceleration_raw_x|acceleration.raw.x", names(acc), value=T), 
#                                axisY=grep("acceleration_raw_y|acceleration.raw.y", names(acc), value=T), 
#                                axisZ=grep("acceleration_raw_z|acceleration.raw.z", names(acc), value=T), 
#                                timeCol="timestamp"){
#   #Continue only if number of acc axes doesn't vary within the same individual
#   if(nrow(acc)==0){stop("The acc dataset has 0 observations.")
#   }else if(nrow(acc)>0){
#     #Split the acc dataframe by burst timestamp, calculate mean vedba per burst and return it as row of the new vedba df
#     burst_ls <- split(acc, as.character(acc[,timeCol]))
#     accDf_vedba <- as.data.frame(rbindlist(lapply(burst_ls, function(burst){
#       n_samples_per_axis <- nrow(burst)
#       acc_burst_duration_s <- 1
#       acc_sampl_freq <- n_samples_per_axis/acc_burst_duration_s
#       vedba <- calculateVedba(accX=burst[,axisX], accY=burst[,axisY], accZ=burst[,axisZ])
#       return(data.frame(timestamp=as.POSIXct(as.character(burst[1,timeCol]), format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
#                         acc_event_id=burst[1,accEventCol], #for the event id, as each acc observation has 1, we only take the first of each burst
#                         n_samples_per_axis=n_samples_per_axis, acc_burst_duration_s=acc_burst_duration_s, acc_sampl_freq=acc_sampl_freq,
#                         meanVedba=mean(vedba, na.rm=T), cumVedba=sum(vedba, na.rm=T),
#                         stringsAsFactors = F))
#     })))
#     return(accDf_vedba)
#   }
# }

