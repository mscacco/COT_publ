
#__________________________
## Species' body mass ####
#__________________________

setwd("...")

#____________________________________________________
## Format dataset with body masses per species ####

# from Hamish Wilman, Jonathan Belmaker, Jennifer Simpson, Carolina de la Rosa, Marcelo M. Rivadeneira, Walter Jetz. 2014. EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. Ecology 95:2027. http://dx.doi.org/10.1890/13-1917.1
# BodyMass-Value, Body mass (g). For Source Dunning08: Geometric mean of average values provided for both sexes (Dunning08). For Source GenAvg: genus average as provided by other sources
# BodyMass-Source, Source of body mass values, Dunning08 or GenAvg PrimScale: inferred from select primary sources with mass and length data, and mass-length relationships parameterized at family level. GenAvg: genus average. Other: see comments.
# BodyMass-SpecLevel, Indicates whether body mass values are based on species-level data, binary 1: based on species level data; 0: inferred from genus or family typical values
# When using the body mass data, please separately cite the respective main sources (two existing digital datasets: Smith et al 2003, Dunning 2007).

# Instead of importing the table we read the lines as there are some rows with different numbers of columns
## IMPORT THE BIRDS DATA
tbl <- readLines("./DataAvailable/BodySizes_perSpecies_dataFrom_WilmanJetz2014/BirdFuncDat.txt")
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
tbl <- readLines("./DataAvailable/BodySizes_perSpecies_dataFrom_WilmanJetz2014/MamFuncDat.txt")
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
acc_studies <- read.csv("./DataAvailable/GpsAcc_MovebankDatasets/ACC_StudiesList_BirdsBats_June2022_perSpecies.csv")

# Check if all species in our list of studies have a corresponding body mass values in the species_info table
table(unique(acc_studies$species) %in% speciesInfo$species_latin)
unique(acc_studies$species[which(!acc_studies$species %in% speciesInfo$species_latin)])

# Some species occurr in our studies list with two different names, so we add both of them to the matching species list
c("Chlamydotis undulata","Chlamydotis macqueenii") %in% acc_studies$species
c("Milvus lineatus","Milvus migrans") %in% acc_studies$species
c("Hieraaetus fasciatus","Aquila fasciata") %in% acc_studies$species

# Build a dataset with the replacement names for those species for which the names don't match
# On the left the names in our study list (matching name)
# On the right the latin name in the downloaded body mass dataset (body mass name)
matchingSpeciesDF <- data.frame(matrix(c(
  "Aquila fasciata","Aquila fasciatus",
  "Hieraaetus fasciatus","Aquila fasciatus",
   "Ardea alba","Casmerodius albus",
  "Torgos tracheliotus","Torgos tracheliotos",
  "Theristicus branickii","Theristicus melanopis",
  "Tringa incana","Heteroscelus incanus",
  "Onychoprion fuscatus","Sterna fuscata",
  "Spizaetus bartelsi","Nisaetus bartelsi",
  "Apus melba","Tachymarptis melba",
  "Anthropoides virgo","Grus virgo",
  "Chlamydotis undulata","Chlamydotis undulata",
  "Chlamydotis macqueenii","Chlamydotis undulata",
  "Milvus migrans","Milvus migrans",
  "Milvus lineatus","Milvus migrans",
  "Hieraaetus wahlbergi","Aquila wahlbergi",
  "Aquila spilogaster","Hieraaetus spilogaster",
  "Larus scoresbii","Leucophaeus scoresbii"
), ncol=2, byrow=T), stringsAsFactors=F)

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

# Add body mass values to the study list to know how many studies we have for different ranges of body masses 
acc_studies_bodyMass <- merge(x=acc_studies, y=speciesInfo[,c("matchingSpeciesName","BodyMass_value","BodyMass_source")], by.x="species", by.y="matchingSpeciesName", all.x=T)
table(is.na(acc_studies_bodyMass$BodyMass_value))
unique(acc_studies_bodyMass$species[which(acc_studies_bodyMass$BodyMass_value<500)])

write.csv(acc_studies_bodyMass, file="./DataAvailable/GpsAcc_MovebankDatasets/allACCstudies_BirdsBats_speciesBodyMass_new.csv", row.names = F)


#_____________________
## RADAR species ####
# Now do the same also for the species in the RADAR studies

speciesInfo_gps <- read.csv("./DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy.csv")
radarSpecies <- read.csv("./DataAvailable/RADARdata_YuvalNir/midjuly_midseptember_Desert_Med_Species-Final.csv")

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
