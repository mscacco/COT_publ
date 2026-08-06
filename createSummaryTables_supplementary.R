
#_____________________________
# Create SUMMARY DATASETS ####
#_____________________________

library(lubridate)

# Set as wd the path where you stored the content downloaded form the Edmond repository
setwd("...")
finalDF <- readRDS("BiologgingData/FinalDf_perPoint_VedbaGs_flappingProbs_ENV.rds")

#______________________________________________________
## Study summary - info per study for supplementary ----

allStudiesLs <- split(finalDF, finalDF$study.name)

summaryInfos_allStudies <- do.call(rbind, lapply(allStudiesLs, function(st){
  return(data.frame(study.name=unique(st$study.name),
                    study.id=unique(st$study.id),
                    deviceType=unique(st$deviceType),
                    species=paste(unique(st$species), collapse="|"),
                    n_indiv=length(unique(st$individual.local.identifier)),
                    n_commutingLocations=nrow(st),
                    n_commutingSegments=length(unique(st$track_flight_id)),
                    accSamplFreq=paste(unique(st$acc_sampl_freq_per_axis), collapse="|"),
                    accBurstDurat_sec=paste(unique(st$acc_burst_duration_s), collapse="|"),
                    med_GPStimeLag_min=median(st$timeLag_min)
  ))
}))
row.names(summaryInfos_allStudies) <- 1:nrow(summaryInfos_allStudies)

write.csv(summaryInfos_allStudies, file="DataFinalSummary/4supplementary_finalSummaryDescription_perStudy_Dec2024.csv", row.names = F)


#_______________________________________
## Species summary for supplementary ----

allSpeciesLs <- split(finalDF, finalDF$species)

speciesDf <- do.call(rbind, lapply(allSpeciesLs, function(sp)try({
  return(data.frame(species=unique(sp$species),
                    Body_mass_kg=unique(sp$Body_mass_kg),
                    deviceType=paste(unique(sp$deviceType), collapse = "|"),
                    totNlocs=nrow(sp),
                    NtrackingDays=length(unique(date(sp$timestamp))),
                    Nstudies=length(unique(sp$study.name)),
                    Nindiv=length(unique(sp$individual.local.identifier)),
                    avg_stepLength_m=mean(sp$stepLength_m, na.rm=T),
                    avg_grSpeed_ms=mean(sp$groundSpeed_ms, na.rm=T),
                    StDev_grSpeed_ms=sd(sp$groundSpeed_ms, na.rm=T),
                    med_timeLag_min=median(sp$timeLag_min, na.rm=T),
                    min_timeLag_min=min(sp$timeLag_min, na.rm=T),
                    tot_distCovered_km=sum(sp$stepLength_m, na.rm=T)/1000,
                    tot_trackDurationAnalysed_h=sum(sp$timeLag_min, na.rm=T)/60,
                    med_accSamplFreq=median(sp$acc_sampl_freq_per_axis),
                    med_accBurstDur_s=median(sp$acc_burst_duration_s),
                    avg_probFlap=mean(sp$flapping_prob, na.rm=T)
  ))
})))
rownames(speciesDf) <- 1:nrow(speciesDf)

# add species names matching with phylogenetic tree
speciesDf$species_phy <- speciesDf$species
speciesDf$species_phy[which(speciesDf$species == "Anthropoides virgo")] <- "Grus virgo"

# save species summary dataset
# This is the final list of species that I will also use to manually associate wing loading values
write.csv(speciesDf, "DataFinalSummary/4supplementary_finalSummaryDescription_perSpecies.csv", row.names = F)

# save species list for subsetting the phylogenetic tree
write.table(as.data.frame(speciesDf$species_phy), file="Phylogeny/speciesList_forPhyloTree.txt", 
            eol="\r\n", quote=F, row.names = F, col.names=F)

