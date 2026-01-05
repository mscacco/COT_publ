
#_____________________________
# Create SUMMARY DATASETS ####
#_____________________________

library(lubridate)

setwd("...")
finalDF <- readRDS("DataFinalSummary/FinalDf_perPoint_VedbaGs_flappingProbs.rds")

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
write.csv(speciesDf, "DataFinalSummary/4supplementary_finalSummaryDescription_perSpecies_Dec2024.csv", row.names = F)

# save species list for subsetting the phylogenetic tree
write.table(as.data.frame(speciesDf$species_phy), file="Phylogeny/speciesList_forPhyloTree.txt", 
            eol="\r\n", quote=F, row.names = F, col.names=F)

#_________________________________________________
## Add wing area and wing loading per species ----

setwd("...")

### Load species list with manually added wing morphology info. Info from JEB2019, BB2001, Ostrich1989 and PLOS2007 were added manually from papers. ----
speciesDf <- read.csv("DataFinalSummary/4supplementary_finalSpeciesList_wingMorphology_Feb2025_final.csv")
head(speciesDf)
# for the wing span range obtained from the birds of the world, take the middle value of the range
speciesDf$wingSpanAverage_birdsOfWorld_cm <- sapply(sapply(strsplit(speciesDf$wingSpan_birdsOfWorld_cm, "-"), as.numeric), mean)

### Load wing measurements (lift generation area and wing area) from Malo & Mata, Ecology and Evolution 2021 ----
mm2021 <- read.table("DataAvailable/MorphologicalMeasurements/WingMorphologyMeasurements/MaloMata2021_WingData.tab", sep="\t", dec=",", na.strings = c("N/A",""))
names(mm2021) <- mm2021[2,]
mm2021 <- mm2021[-c(1:2),c("Scientific name","Lift generation area (in m2)","Wing area (S, in m2)","Source","Wing span (cm)","Source wing span", "Number of samples")]
head(mm2021)

mm2021$`Scientific name` <- gsub("\t|=|-","", mm2021$`Scientific name`)
mm2021$`Lift generation area (in m2)` <- as.numeric(gsub(",",".", mm2021$`Lift generation area (in m2)`))
mm2021$`Wing area (S, in m2)` <- as.numeric(gsub(",",".", mm2021$`Wing area (S, in m2)`))
mm2021$`Wing span (cm)` <- as.numeric(gsub(",",".", mm2021$`Wing span (cm)`))

table(speciesDf$species %in% mm2021$`Scientific name`)
length(grep(paste(speciesDf$species, collapse="|"), mm2021$`Scientific name`, value=T))

speciesDf[,c("WingSpan_cm_MM2021", "WingSpan_MM2021_origSource", "LiftGenArea_cm2_MM2021", "WingArea_cm2_MM2021", "WingArea_MM2021_origSource", "Nsamples_MM2021")] <- NA
for(i in 1:nrow(speciesDf)){
  j <- grep(speciesDf$species[i], mm2021$`Scientific name`)
  if(length(j)>0){
    speciesDf[i, "LiftGenArea_cm2_MM2021"] <- mm2021[j, "Lift generation area (in m2)"]*10000
    speciesDf[i, "WingArea_cm2_MM2021"] <- mm2021[j, "Wing area (S, in m2)"]*10000
    speciesDf[i, "WingSpan_cm_MM2021"] <- mm2021[j, "Wing span (cm)"]
    speciesDf[i, c("WingArea_MM2021_origSource","WingSpan_MM2021_origSource","Nsamples_MM2021")] <- mm2021[j, c("Source","Source wing span","Number of samples")]
  }}
head(speciesDf)
summary(speciesDf)

# check the wing span value side by side, they are very similar for species for which we have values in both
speciesDf[,c("wingSpanAverage_birdsOfWorld_cm","WingSpan_cm_MM2021")]
# In order to calculate the wing area proxy, let's keep the values of MM and add the average values from birds of the world where MM are missing
speciesDf$WS_final_cm <- speciesDf$WingSpan_cm_MM2021
speciesDf$WS_final_cm[is.na(speciesDf$WingSpan_cm_MM2021)] <- speciesDf$wingSpanAverage_birdsOfWorld_cm[is.na(speciesDf$WingSpan_cm_MM2021)]
speciesDf[,c("species","species_phy","wingSpanAverage_birdsOfWorld_cm","WingSpan_cm_MM2021","WS_final_cm")]

### Load wing measurements from AVONET ----
avonet <- read.csv("DataAvailable/MorphologicalMeasurements/WingMorphologyMeasurements/AVONET_Supp1_Avonet3BirdTree.csv")
table(speciesDf$species_phy %in% avonet$Species3)
speciesDf$species_phy[!speciesDf$species_phy %in% avonet$Species3]

speciesDf <- merge(speciesDf, avonet[,c("Species3","Wing.Length","Kipps.Distance","Secondary1","Hand.Wing.Index")], by.x="species_phy", by.y="Species3", all.x=T)
speciesDf <- speciesDf %>% rename(WingLength_mm=Wing.Length) %>% rename(Secondary1_mm=Secondary1)  %>% rename(KippsDistance_mm=Kipps.Distance)
head(speciesDf)

### Apply ellipse folded-wing area estimation from Hellen 2023 Ecology and Evolution "New methods for estimating the total wing area of birds" ----
ellipse <- function(WL, S1, WS){(1/2 * pi * WL * S1) + ((WS - (2 * WL)) * S1)}
speciesDf$wingArea_ellipse_cm2 <- ellipse(WL=speciesDf$WingLength_mm/10,
                                      S1=speciesDf$Secondary1_mm/10, 
                                      WS=speciesDf$WS_final_cm)
head(speciesDf)

### Calculate wing loading
speciesDf$wingLoading_kgm2 <-  speciesDf$Body_mass_kg / (speciesDf$wingArea_ellipse_cm2 * 1e-4)
range(speciesDf$wingLoading_kgm2, na.rm=T)
speciesDf[,c("species", "wingLoading_kgm2")]

write.csv(speciesDf, "DataFinalSummary/4supplementary_finalSpeciesList_wingMorphology_Feb2025_final.csv", row.names = F)

# Final notes: data from ICB2022 are the most different, calculated in a different way.
# All others are quite similar for those species for which we have measurements from more than one source. 
# JEB2019 for pelican quite smaller than the other ones with sample size 1. 
# BB2001, Ostrich1989, PLOS2007 and MM2021 (wing area) seem very compatible and could be maybe averaged for those species for which multiple measurements are available?


