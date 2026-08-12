
#__________________________
## Species' morphology ####
#__________________________

# In order to run this script, you will need:
# - the wind morphology dataset compiled by the authors, and available in the Edmond repository "InputData/wingMorphology_perSpecies_Feb2025.csv"
# - wing measurements (lift generation area and wing area) available from Malo & Mata, Ecology and Evolution 2021
# - the Avonet dataset, available from "AVONET: morphological, ecological and geographical data for all birds". Ecol. Lett. 25, 581–597 (2022)
# - the body mass dataset, available from "EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals". Ecology 95:2027 (2014) 

# Set as wd the path where you stored the content downloaded from the Edmond repository
setwd("...")

#____________________________________________
## 1. Compile wing morphology per species ####

# ! Note that the file "InputData/wingMorphology_perSpecies_Feb2025.csv" contains both the input and the output of this section
# as the columns output of this script get simply added to the same file

# Store the path where you downloaded all additional wing datasets
wm_path <- "..."


### Load species list with manually added wing morphology info. Info from JEB2019, BB2001, Ostrich1989 and PLOS2007 were added manually from papers. ----
speciesDf <- read.csv("InputData/wingMorphology_perSpecies_Feb2025.csv")
head(speciesDf)
# for the wing span range obtained from the birds of the world, take the middle value of the range
speciesDf$wingSpanAverage_birdsOfWorld_cm <- sapply(sapply(strsplit(speciesDf$wingSpan_birdsOfWorld_cm, "-"), as.numeric), mean)

### Load wing measurements (lift generation area and wing area) from Malo & Mata, Ecology and Evolution 2021 ----
mm2021 <- read.table(paste0(wm_path,"/MaloMata2021_WingData.tab"), sep="\t", dec=",", na.strings = c("N/A",""))
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
avonet <- read.csv(paste0(wm_path,"/AVONET_Supp1_Avonet3BirdTree.csv"))
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

write.csv(speciesDf, "InputData/wingMorphology_perSpecies_Feb2025.csv", row.names = F)


#____________________________________________________
## 2. Format dataset with body masses per species ####

# !! Note: in order to run this script, the body mass dataset needs to be downloaded separate from:
# http://dx.doi.org/10.1890/13-1917.1 
# Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, Marcelo M. Rivadeneira, Walter Jetz. 2014. 
# "EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals". Ecology 95:2027.
# and stored in a path of choice:

bm_path <- "..."

# BodyMass-Value, Body mass (g). For Source Dunning08: Geometric mean of average values provided for both sexes (Dunning08). For Source GenAvg: genus average as provided by other sources
# BodyMass-Source, Source of body mass values, Dunning08 or GenAvg PrimScale: inferred from select primary sources with mass and length data, and mass-length relationships parameterized at family level. GenAvg: genus average. Other: see comments.
# BodyMass-SpecLevel, Indicates whether body mass values are based on species-level data, binary 1: based on species level data; 0: inferred from genus or family typical values
# When using the body mass data, please separately cite the respective main sources (two existing digital datasets: Smith et al 2003, Dunning 2007).

# Instead of importing the table we read the lines as there are some rows with different numbers of columns
## IMPORT THE BIRDS DATA
tbl <- readLines(paste0(bm_path, "BirdFuncDat.txt"))
colNames <- unlist(strsplit(tbl[1], "\t")) #isolate headers
tbl <- tbl[-1] #Remove them from the rest of the table
# Check that all entries have the same number of columns (separator is tab \t)
table(sapply(strsplit(tbl, "\t"), length)) # all entries except a few are missing the "Record-Comment" column
# Remove entries with < 39 cols and add an extra empty column to all entries which miss that
cols39 <- sapply(strsplit(tbl, "\t"), length)==39
cols40 <- sapply(strsplit(tbl, "\t"), length)==40
ls_allCols <- strsplit(tbl, "\t")[which(cols40)]
ls_missCols <- strsplit(tbl, "\t")[which(cols39)]
ls_missCols <- lapply(ls_missCols, function(x){
  return(c(x, ""))
})
# Now that they have the same number of columns bind the two lists and rbind them into a dataframe, adding the header
birdsInfo <- do.call(rbind, c(ls_missCols, ls_allCols))
colnames(birdsInfo) <- colNames
birdsInfo <- as.data.frame(birdsInfo, stringsAsFactors=F)
str(birdsInfo)
# Subset only columns of interest and save table
birdsInfo <- birdsInfo[,c("BLFamilyLatin","Scientific","English","BodyMass-Value","BodyMass-Source","BodyMass-SpecLevel")]
names(birdsInfo) <- c("FamilyLatin","species_latin","species_english","BodyMass_value","BodyMass_source","BodyMass_specLevel")

write.csv(birdsInfo, "./DataAvailable/birdSpecies_bodyMasses.csv", row.names = F)


#_________________________________
## IMPORT THE BATS (MAMMALS) DATA
tbl <- readLines(paste0(bm_path, "MamFuncDat.txt"))
colNames <- unlist(strsplit(tbl[1], "\t")) #isolate headers
tbl <- tbl[-1] #Remove them from the rest of the table
# Check that all entries have the same number of columns (separator is tab \t)
table(sapply(strsplit(tbl, "\t"), length)) # all entries except a few are missing the "Record-Comment" column
# Remove entries with < 26 cols and add an extra empty column to all entries which miss that
cols25 <- sapply(strsplit(tbl, "\t"), length)==25
cols26 <- sapply(strsplit(tbl, "\t"), length)==26
ls_allCols <- strsplit(tbl, "\t")[which(cols26)]
ls_missCols <- strsplit(tbl, "\t")[which(cols25)]
ls_missCols <- lapply(ls_missCols, function(x){
  return(c(x, ""))
})
# Now that they have the same number of columns bind the two lists and rbind them into a dataframe, adding the header
mammInfo <- do.call(rbind, c(ls_missCols, ls_allCols))
colnames(mammInfo) <- colNames
mammInfo <- as.data.frame(mammInfo, stringsAsFactors=F)
str(mammInfo)
# Subset only columns of interest and save table
mammInfo <- mammInfo[!(is.na(mammInfo$Scientific) | mammInfo$Scientific==""),]
mammInfo <- mammInfo[,c("MSWFamilyLatin","Scientific","BodyMass-Value","BodyMass-Source","BodyMass-SpecLevel")]
names(mammInfo) <- c("FamilyLatin","species_latin","BodyMass_value","BodyMass_source","BodyMass_specLevel")
mammInfo$species_english <- NA

# save table
write.csv(mammInfo, "./DataAvailable/mammSpecies_bodyMasses.csv", row.names = F)


## Rbind the two datasets with the species info of birds and bats
speciesInfo <- rbind(birdsInfo, mammInfo)

# Remove rows with empty values from body mass or rows with non-numeric body mass values
anyNA(speciesInfo$species_latin)
anyNA(speciesInfo$BodyMass_value)
speciesInfo$species_latin[which(is.na(speciesInfo$BodyMass_value))]
speciesInfo <- speciesInfo[which(!speciesInfo$species_latin==""),]
grep("[a-zA-Z]",speciesInfo$BodyMass_value)

# Set body mass to numeric
class(speciesInfo$BodyMass_value)
speciesInfo$BodyMass_value <- as.numeric(speciesInfo$BodyMass_value)
str(speciesInfo)

# Save table
write.csv(speciesInfo,  file="./DataAvailable/BodyMassInfos_allBirdBatsSpecies.csv", row.names = F)


#_______________________
## GPS-ACC species ####
## Import list of Movebank studies to download to check that the taxonomy names match for later merging

speciesInfo <- read.csv("./DataAvailable/BodyMassInfos_allBirdBatsSpecies.csv")
acc_studies <- read.csv("./DataAvailable/ACC_StudiesList_BirdsBats_June2022_perSpecies.csv")

# Check if all species in our list of studies have a corresponding body mass values in the species_info table
table(unique(acc_studies$species) %in% speciesInfo$species_latin)
unique(acc_studies$species[which(!acc_studies$species %in% speciesInfo$species_latin)])

# Build a dataset with the replacement names for those species for which the names don't match
# On the left the names in our study list (matching name)
# On the right the latin name in the downloaded body mass dataset (body mass name)
matchingSpeciesDF <- data.frame(matrix(
  c("Anthropoides virgo","Grus virgo"), ncol=2, byrow=T), stringsAsFactors=F)
names(matchingSpeciesDF) <- c("matchingSpeciesName","nameInBodyMassDF")

# check that the two names match the respective datasets
table(matchingSpeciesDF$matchingSpeciesName %in% acc_studies$species)
table(matchingSpeciesDF$nameInBodyMassDF %in% speciesInfo$species_latin)
# add the new matching name to the body mass DF (speciesInfo dataset)
speciesInfo <- merge(x=speciesInfo, y=matchingSpeciesDF, by.x="species_latin", by.y="nameInBodyMassDF", all.x=T)
# Fill the NAs in the matching column (corresponding to the names that didn't change) with the original name
speciesInfo$matchingSpeciesName[is.na(speciesInfo$matchingSpeciesName)] <- speciesInfo$species_latin[is.na(speciesInfo$matchingSpeciesName)]
# Check that those that are different are only those that were missing
table(speciesInfo$species_latin %in% speciesInfo$matchingSpeciesName)

# Save it
write.csv(speciesInfo,  file="./DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv", row.names = F)

#_____________________
## RADAR species ####
# Now do the same also for the species in the RADAR studies

speciesInfo_gps <- read.csv("./DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv")
radarSpecies <- read.csv("./RadarData/SpeciesList_midjuly_midseptember_Desert_Med_Species-Final.csv")

radarSpecies$species_latin <- paste(radarSpecies$Genus, speciesInfos$species, sep=" ")
length(unique(radarSpecies$species_latin))
unique(df$July)

table(radarSpecies$species_latin %in% speciesInfo_gps$matchingSpeciesName)
radarSpecies$species_latin[!radarSpecies$species_latin %in% speciesInfo_gps$matchingSpeciesName]

matchingSpeciesDF <- data.frame(matrix(c(
  "Carpospiza brachydactyla","Petronia brachydactyla",
  "Cercotrichas galactotes","Erythropygia galactotes",
  "Curruca communis","Sylvia communis",
  "Curruca crassirostris","Sylvia hortensis", # just recently split in two species, one is western one is easter orphean warbler
  "Curruca curruca","Sylvia curruca",
  "Curruca nisoria","Sylvia nisoria",
  "Curruca ruppeli","Sylvia rueppelli",
  "Iduna pallida","Hippolais pallida",
  "Oenanthe melanoleuca","Oenanthe hispanica" # by some considered conspecific (one is easter oen is western black-eared wheatear)
), ncol=2, byrow=T), stringsAsFactors=F)

names(matchingSpeciesDF) <- c("matchingSpeciesName_radar","nameInBodyMassDF")

# check that the two names match the respective datasets
table(matchingSpeciesDF$matchingSpeciesName_radar %in% radarSpecies$species_latin)
table(matchingSpeciesDF$nameInBodyMassDF %in% speciesInfo_gps$matchingSpeciesName)
# add the new matching name to the body mass DF (speciesInfo dataset)
speciesInfo_gps_radar <- merge(x=speciesInfo_gps, y=matchingSpeciesDF, by.x="matchingSpeciesName", by.y="nameInBodyMassDF", all.x=T)
# Fill the NAs in the matching column (corresponding to the names that didn't change) with the original name
head(speciesInfo_gps_radar)
table(is.na(speciesInfo_gps_radar$matchingSpeciesName_radar))
speciesInfo_gps_radar$matchingSpeciesName_radar[is.na(speciesInfo_gps_radar$matchingSpeciesName_radar)] <- speciesInfo_gps_radar$matchingSpeciesName[is.na(speciesInfo_gps_radar$matchingSpeciesName_radar)]

# Check that all species from both GPS and RADAR studies match the names
table(radarSpecies$species_latin %in% speciesInfo_gps_radar$matchingSpeciesName_radar)
table(acc_studies$species %in% speciesInfo_gps_radar$matchingSpeciesName_radar)

# Change column names
names(speciesInfo_gps_radar)[names(speciesInfo_gps_radar)=="matchingSpeciesName_radar"] <- "matchingSpeciesName_gps_radar"
speciesInfo_gps_radar$matchingSpeciesName <- NULL # remove the old matching col
head(speciesInfo_gps_radar)
table(is.na(speciesInfo_gps_radar$matchingSpeciesName_gps_radar))

# Save it
write.csv(speciesInfo_gps_radar,  file="./DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy_gps+radar.csv", row.names = F)
