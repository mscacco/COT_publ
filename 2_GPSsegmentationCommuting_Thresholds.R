
#__________________________________________
# Classification of commuting flights ####
#__________________________________________

#____________________________________________________________
# Explore possible segmentation for commuting segments ####

setwd("...")
dir.create("Plots/SPEEDhistograms")

fls <- list.files("DataProcessed", "dfGpsAcc_geom_10min.RData", full.names=T)

pdf(file="Plots/SPEEDhistograms/speeds_allStudies.pdf")
err <- lapply(fls, function(f){
  studySpeciesId <- paste(strsplit(f, "_")[[1]][2:3], collapse="_")
  load(f) #data.frame object df_geom
  print(paste0(unique(df_geom$study.name)," - ",studySpeciesId))
  df <- df_geom[df_geom$timeLag_min < 70,]
  if(nrow(df)>10){ #Keep the individual only if it has at least 10 observations with timelag < 70
    # Keep only obs with no missing values in these variables
    df <- df[complete.cases(df[,c("groundSpeed_ms","turnAngle")]),]
    if(nrow(df[df$groundSpeed_ms>2,])>10){ #Keep the individual only if it has at least 10 observations above 2m/s
      # myBW <- bw.nrd(abs(df$groundSpeed_ms))
      # return(nmodes(abs(df$groundSpeed_ms), bw=myBW))
      return(hist(df$groundSpeed_ms, breaks="FD", main=studySpeciesId, ylim=c(0,1000), xlim=c(0,25)))
    }
  }
})
dev.off()

pdf(file="Plots/SPEEDhistograms/turningAngles_highSpeeds_allStudies.pdf")
err <- lapply(fls, function(f){
  studySpeciesId <- paste(strsplit(f, "_")[[1]][2:3], collapse="_")
  load(f) #data.frame object df_geom
  print(paste0(unique(df_geom$study.name)," - ",studySpeciesId))
  df <- df_geom[df_geom$timeLag_min < 70,]
  if(nrow(df)>10){ #Keep the individual only if it has at least 10 observations with timelag < 70
    # Keep only obs with no missing values in these variables
    df <- df[complete.cases(df[,c("groundSpeed_ms","turnAngle")]),]
    df_sub <- df[df$groundSpeed_ms > 2.5,]
    if(nrow(df_sub)>10){ #Keep the individual only if it has at least 10 observations above 2m/s
      return(hist(abs(df_sub$turnAngle), breaks="FD", main=studySpeciesId, ylim=c(0,5000), xlim=c(0,180)))
    }
  }
})
dev.off()


pdf(file="Plots/SPEEDhistograms/meanVedba_highSpeed_lowTurnAngle_allStudies.pdf")
err <- lapply(fls, function(f){
  studySpeciesId <- paste(strsplit(f, "_")[[1]][2:3], collapse="_")
  load(f) #data.frame object df_geom
  print(paste0(unique(df_geom$study.name)," - ",studySpeciesId))
  df <- df_geom[df_geom$timeLag_min < 70,]
  if(nrow(df)>10){ #Keep the individual only if it has at least 10 observations with timelag < 70
    # Keep only obs with no missing values in these variables
    df <- df[complete.cases(df[,c("groundSpeed_ms","turnAngle")]),]
    df_sub <- df[df$groundSpeed_ms > 2.5 & df$turnAngle < 100,]
    if(nrow(df_sub)>10){ #Keep the individual only if it has at least 10 observations above 2m/s
      return(hist(df_sub$meanVedba, breaks="FD", main=studySpeciesId, ylim=c(0,5000)))
    }
  }
})
dev.off()


#_____________________________________________________
# Extract commuting segments from complete tracks ####

library(ggplot2)
theme_set(theme_bw())
library(sf)
library(plyr)
library(doParallel)
doParallel::registerDoParallel(2)
library(rnaturalearth)
world_map <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')

setwd("...")

fls <- list.files("DataProcessed", "dfGpsAcc_geom_10min.RData", full.names=T)
dir.create("Plots/finalFlightSegments_thresholdClassification") # to store plots

err <- lapply(fls, function(f)try({
  studySpeciesId <- paste(strsplit(f, "_")[[1]][2:3], collapse="_")
  deviceType <- strsplit(f, "_")[[1]][4]
  load(f) #data.frame object df_geom
  if(any(class(df_geom)=="data.table")){df_geom <- as.data.frame(df_geom)}# make sure is not data.table
  print(paste0(unique(df_geom$study.name)," - ",studySpeciesId))
  # exclude long timelags
  df_geom <- df_geom[which(df_geom$timeLag_min < 70),]

  if(nrow(df_geom) > 10){
    if(nrow(df_geom[which(df_geom$groundSpeed_ms > 1.5),]) > 1){ #Keep the individual only if it has at least 10 observations above 2m/s
      # Split each species df by track id
      indivID <- if("individual.tag.id" %in% names(df_geom)){"individual.tag.id"}else{"individual.local.identifier"}
      indLs <- split(df_geom, df_geom[,indivID])
      
      df_allSegm <- do.call(rbind, lapply(indLs, function(t){
        t <- t[order(t$timestamp),]
        t$flightClass <- "nonFlying"
        t$flightClass[which(abs(t$turnAngle) < 100 & t$groundSpeed_ms > 1.5)] <- "Flying"
        # Assign segment ID to consecutive observations in flight
        flightNum <- rep(1, nrow(t)) # "nonFlying" and NAs get a 1 (cumsum does not work if there are NAs)
        flightNum[which(t$flightClass=="Flying")] <- 0 # only "flying" gets a 0
        t$track_flight_id <- paste0(unique(t[,indivID]),"_",cumsum(flightNum)) # Assign unique ID to the flight segment
        t$track_flight_id[t$flightClass!="Flying"|is.na(t$flightClass)] <- NA # set to NA the segment id that are not classified as flight (either nonFlying or NA)
        return(t)
      }))
      # Keep among the flight segments only those that are minimum 3 observations long (> 10 min if minimum sampling frequency is 5 min, 3 hours with hourly data) 
      longSegments <- names(which(table(df_allSegm$track_flight_id) >= 3))
      df_allSegm$track_flight_id[which(!df_allSegm$track_flight_id %in% longSegments)] <- NA # set to NA the ID of the segments that do not satisfy this requirement
      
      # Save dataset with all segments (full trajectory) but with ID associated to flight segments > 10 observations
      save(df_allSegm, file=paste0("DataProcessed/studyId_",studySpeciesId,"_",deviceType,"_dfGpsAcc_allSegmentsID_thresholdClass.RData"))
      
      # If there are some long enough segments in the dataset:
      # Plot the traj with highlighted final segments
      colFact <- factor(as.numeric(!is.na(df_allSegm$track_flight_id)), levels=c("1","0"))
      if(length(longSegments)>0){
         colLabs <- c("Commuting","Other"); colVals <- c("deeppink","darkgrey")}else{ ##006ad1
          colLabs <- "Other"; colVals <- "darkgrey"; warning("No observations were classified as commuting flight.")
         }
       ggplot() +
        geom_sf(data = world_map, fill="black", colour = "white") +#, extend=F) +
        ylim(c(min(df_allSegm$coords.x2), max(df_allSegm$coords.x2))) +
        xlim(c(min(df_allSegm$coords.x1), max(df_allSegm$coords.x1))) +
        geom_path(df_allSegm, mapping = aes(x=coords.x1, y=coords.x2, group=df_allSegm[,indivID]), col="darkgrey") +
        geom_point(df_allSegm, mapping = aes(x=coords.x1, y=coords.x2, color=colFact, group=df_allSegm[,indivID]), alpha=0.6) +
        scale_color_manual(labels=colLabs, values=colVals, na.value="black") +
        #geom_point(df_allSegm[!is.na(df_allSegm$track_flight_id),], mapping = aes(x=coords.x1, y=coords.x2), alpha=0.6, color="#006ad1") +
        geom_point(df_allSegm[!is.na(df_allSegm$track_flight_id),], mapping = aes(x=coords.x1, y=coords.x2), alpha=0.6, color="deeppink") +
        xlab("") + ylab("") + ggtitle(studySpeciesId) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(col=guide_legend("Flight"))
      ggsave(paste0("Plots/finalFlightSegments_thresholdClassification/selectedCommutingSegments_GPS_blackBG/study_",studySpeciesId,"_",deviceType,"_finalSelectedSegments_threshold_speed15-turn100.png"),
             width=8,height=6,units="in",dpi=300)
    }else if(nrow(df_geom[df_geom$groundSpeed_ms>2,])<=1){stop("Less than 1 observations with ground speed > 1.5 m/s, study will be ignored.")}
  }else if(nrow(df_geom)<=10){stop("Less than 10 observations with timelags < 70 minutes, study will be ignored.")}
}))

# Check errors
is.error <- function(x) inherits(x, "try-error")
errors <- vapply(err, is.error, logical(1))
err[vapply(err, is.error, logical(1))]
fls[errors]

