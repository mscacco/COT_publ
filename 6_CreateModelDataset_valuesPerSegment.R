
#_________________________________________________
# CREATE SEGMENT DATASET FOR MODELS (step 6) ####
#_________________________________________________

library(lubridate)
library(data.table)

#setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

# Load the final dataset (point by point) with all the annotated atmospheric variables from step 4(A-C)
finalDF <- readRDS("DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV_Feb2025.rds")

#__________________________________________________
# Calculate a summary information per SEGMENT ####

# Segments have a variable duration, mostly between 12 min and 10 hours (so segment duration will have to be in the model)
# Ndays_segment <- aggregate(timestamp~track_flight_id, data=finalDF, FUN=function(x)length(unique(date(x))))
# summary(Ndays_segment[,2][Ndays_segment[,2] > 1]) # between 1 and 6 days
Nmin_segment <- aggregate(timeLag_min~track_flight_id, data=finalDF, FUN=sum)
summary(Nmin_segment[,2]/60) #segment duration in hours
hist(Nmin_segment[,2]/60)

segmentLs <- split(finalDF, finalDF$track_flight_id)
segmentDf <- as.data.frame(rbindlist(lapply(segmentLs, function(sp){
  sp <- sp[order(sp$timestamp),]
  return(data.frame(segmentID=unique(sp$track_flight_id),
                    individualID=unique(sp$individual.local.identifier),
                    studyName=unique(sp$study.name),
                    studyID=unique(sp$study.id),
                    species=unique(sp$species),
                    Body_mass_g=unique(sp$BodyMass_value),
                    Body_mass_kg=unique(sp$Body_mass_kg),
                    deviceType=unique(sp$deviceType),
                    dataSource=unique(sp$dataSource),
                    totNlocs=nrow(sp),
                    segm_timeStart=sp$timestamp[1],
                    segm_timeEnd=sp$timestamp[nrow(sp)],
                    segm_timeDuration_h=as.numeric(difftime(sp$timestamp[nrow(sp)], sp$timestamp[1], units="hour")),
                    avg_stepLength_m=mean(sp$stepLength_m, na.rm=T),
                    avg_grSpeed_ms=mean(sp$groundSpeed_ms, na.rm=T),
                    avg_timeLag_min=mean(sp$timeLag_min, na.rm=T),
                    min_timeLag_min=min(sp$timeLag_min, na.rm=T),
                    StDev_grSpeed_ms=sd(sp$groundSpeed_ms, na.rm=T),
                    tot_distCovered_km=sum(sp$stepLength_m, na.rm=T)/1000,
                    tot_trackDurationAnalysed_h=sum(sp$timeLag_min, na.rm=T)/60,
                    # env metrics
                    avg_dem=mean(sp$dem, na.rm=T),
                    avg_oroUplift=mean(sp$Orographic_uplift_potential, na.rm=T),
                    avg_thermUplift=mean(sp$`w_star_Thermal_uplift_potential_(W*)`, na.rm=T),
                    avg_heatFlux=mean(sp$ishf_Instantaneous_surface_sensible_heat_flux, na.rm=T),
                    avg_temp_1000m=mean(sp$t_Temperature, na.rm=T),
                    avg_windSpeed_1000m=mean(sp$windSpeed_1000m, na.rm=T),
                    avg_windSupp_1000m=mean(sp$windSupport_1000m, na.rm=T),
                    avg_crossWind_1000m=mean(sp$crossWind_1000m, na.rm=T),
                    avg_airspeed_1000m=mean(sp$airspeed_1000m, na.rm=T),
                    # acc metrics
                    avg_VedbaGs=mean(sp$meanVedba_Gs, na.rm=T),
                    med_accSamplFreq=median(sp$acc_sampl_freq_per_axis, na.rm=T),
                    med_accBurstDur_s=median(sp$acc_burst_duration_s, na.rm=T),
                    # active flight metrics
                    avg_probFlap=mean(sp$flapping_prob, na.rm=T),
                    # passive metrics
                    avg_probPass=mean((1-sp$flapping_prob), na.rm=T),
                    # proportion of flapping per segment, as average of 0s and 1s
                    avg_propFlap_multimode=mean(sp$flapping_bin, na.rm=T))
  )
})))
# Round the numeric variables
colsToRound <- names(which(sapply(segmentDf, is.numeric)))
colsToRound <- colsToRound[!colsToRound %in% "totNlocs"]
segmentDf[,colsToRound] <- round(segmentDf[,colsToRound], 3)
# prop pass and prop flap sum to 1
summary(segmentDf$avg_probFlap + segmentDf$avg_probPass)
# some summary stats
summary(segmentDf$segm_timeDuration_h)
q <- quantile(segmentDf$segm_timeDuration_h, seq(0.9,1,0.0001))
length(unique(segmentDf$species[segmentDf$segm_timeDuration_h <=10]))
summary(segmentDf$avg_timeLag_min)
summary(segmentDf$totNlocs)

# Filter segments longer than 10 hours
segmentDf_sub <- segmentDf[segmentDf$segm_timeDuration_h <= 10,]
saveRDS(segmentDf_sub, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Feb2025_max10hours.rds")

# Summarise some information fo the paper
summary(segmentDf_sub$segm_timeDuration_h)
summary(segmentDf_sub$avg_timeLag_min)
summary(segmentDf_sub$totNlocs)

# Number of studies and species in Movebank
length(unique(allSegmDfs$studyName[allSegmDfs$dataSource=="Movebank"]))
length(unique(allSegmDfs$species[allSegmDfs$dataSource=="Movebank"]))
length(unique(allSegmDfs$individualID[allSegmDfs$dataSource=="Movebank"]))

# Number of studies and species from external sources
length(unique(allSegmDfs$studyName[allSegmDfs$dataSource!="Movebank"]))
length(unique(allSegmDfs$species[allSegmDfs$dataSource!="Movebank"]))
length(unique(allSegmDfs$individualID[allSegmDfs$dataSource!="Movebank"]))

# Devices
table(allSegmDfs$deviceType)

#__________________________________________________________________________________
# Check that species names match names in Phylogenetic Tree from birdtree.org ####

treeDf <- read.csv("Phylogeny/2012-03-04206D-master_taxonomy.csv")

# check if names match with those in the phylogenetic tree
unique(segmentDf_sub$species)[!unique(segmentDf_sub$species) %in% treeDf$Scientific]
# these are the corresponding names for those species
"Grus virgo" %in% treeDf$Scientific

# we add a column with matching names
segmentDf_sub$species_phy <- segmentDf_sub$species
segmentDf_sub$species_phy[which(segmentDf_sub$species == "Anthropoides virgo")] <- "Grus virgo"

# now they all match (except the 2 bat species)
unique(segmentDf_sub$species_phy)[!unique(segmentDf_sub$species_phy) %in% treeDf$Scientific]

# re-save dataset to use in step 6B
#write.csv(segmentDf_sub, "DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Dec2024_max10hours.csv", row.names = F)
saveRDS(segmentDf_sub, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Feb2025_max10hours.rds")

