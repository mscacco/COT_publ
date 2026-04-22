
#_______________________________________
# CREATE SEGMENT DATASET FOR MODELS ####
#_______________________________________

library(lubridate)
library(data.table)

setwd("...")

# Load the final dataset (point by point) with all the annotated atmospheric variables from step 5(A-B)
finalDF <- readRDS("DataFinalSummary/FinalDf_perPoint_VedbaGs_flappingProbs_ENV.rds")

#__________________________________________________
# Calculate a summary information per SEGMENT ####

# Segments have a variable duration, mostly between 12 min and 10 hours (so segment duration will have to be in the model)
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
saveRDS(segmentDf_sub, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Feb2025_max10hours.rds")


#___________________
# MAP for FIG 1 ####
#___________________

#________________________________________________________________________________________________________
# Use the point dataset and the segment to map the trajectories of all individuals available for analysis

library(rnaturalearth)
library(ggplot2)
library(scales)
library(colorspace)
library(viridisLite)
library(sf)
library(terra)
library(dplyr)
library(svglite)
# extent to bbox enlarging function:
ext_to_bb <- function(x, f){ # f being a multiplier of how much, proportionally, we want to extent bb
  e <- ext(x)
  dx <- (xmax(e) - xmin(e)) * f; dy <- (ymax(e) - ymin(e)) * f
  e <- ext(xmin(e) - dx, xmax(e) + dx, ymin(e) - dy, ymax(e) + dy)
  bb <- st_as_sfc(st_bbox(e))
  st_crs(bb) <- st_crs(x)
  return(bb)
}

# filter df with all GPS fixes to contain only the trajectories of individuals selected for analyses
pointDf <- readRDS("DataFinalSummary/FinalDf_perPoint_VedbaGs_flappingProbs_ENV.rds")
segmentDf <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Feb2025_max10hours.rds")
pointDf_sub <- pointDf[pointDf$individual.local.identifier %in% unique(segmentDf$individualID),]
pointDf_sub <- pointDf_sub[pointDf_sub$track_flight_id %in% unique(segmentDf$segmentID),]

# transform into sf object
points_sf <- st_as_sf(pointDf_sub, coords = c("location.long", "location.lat"), crs = 4326) 

# Thin all trajectories and turn into lines for faster plotting
points_sf <- arrange(points_sf, individual.local.identifier, tag.local.identifier, timestamp)
points_sf_thin <- group_by(points_sf, individual.local.identifier, tag.local.identifier) %>% slice(seq(1, n(), by = 4))

lines_sf <- points_sf_thin %>%
  arrange(individual.local.identifier, tag.local.identifier, track_flight_id, timestamp) %>%   # Right sorting!!
  group_by(individual.local.identifier, tag.local.identifier, track_flight_id, study.name, species) %>%
  summarise(geometry = st_sfc(
    st_linestring(as.matrix(st_coordinates(geometry))),
    crs = st_crs(points_sf)),
    .groups = "drop") %>%
  st_as_sf()
lines_sf <- lines_sf %>% mutate(study_species = paste(study.name, species, sep = "_"))

length(unique(segmentDf$segmentID)) == nrow(lines_sf) # 850 individuals
length(unique(lines_sf$individual.local.identifier)) # 850 individuals
length(unique(lines_sf$study.name)) # 64 studies
length(unique(lines_sf$study_species)) # 72 combinations study - species

# Select trajectories for zoomed in panel in the Americas at full resolution (not thinned)
westCoast <- points_sf[st_coordinates(points_sf)[,1] < -70,]
plot(st_geometry(westCoast))

lines_inset <- westCoast %>%
  arrange(individual.local.identifier, tag.local.identifier, timestamp) %>%   # Right sorting!!
  group_by(individual.local.identifier, tag.local.identifier, study.name, species) %>%
  summarise(geometry = st_sfc(st_linestring(as.matrix(st_coordinates(geometry))), crs = st_crs(points_sf)),
            .groups = "drop") %>%
  st_as_sf()
lines_inset <- lines_inset %>% mutate(study_species = paste(study.name, species, sep = "_"))


westCoast_comm <- westCoast[westCoast$track_flight_id %in% unique(segmentDf$segmentID),]
commuting_inset <- westCoast_comm %>%
  arrange(individual.local.identifier, tag.local.identifier, track_flight_id,  timestamp) %>%   # Right sorting!!
  group_by(individual.local.identifier, tag.local.identifier, track_flight_id, study.name, species) %>%
  summarise(geometry = st_sfc(st_linestring(as.matrix(st_coordinates(geometry))), crs = st_crs(points_sf)),
            .groups = "drop") %>%
  st_as_sf()
commuting_inset <- commuting_inset %>% mutate(study_species = paste(study.name, species, sep = "_"))

#_____________
# Background map and colors

dir.create("Revision/NewSupplFigures/globalMap")

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$continent != "Antarctica", ]

n <- length(unique(lines_sf$study_species))
cols <- scales::hue_pal(h = c(0, 360), c = 60, l = 65)(n)
hcl <- as(hex2RGB(cols), "polarLUV")@coords
# remove magenta tone (roughly hue 300°)
keep <- !(hcl[,3] > 310 | hcl[,3] < 15)
hcl <- hcl[keep, ]
cols_clean <- hex(polarLUV(
  L = pmin(95, hcl[,1] + 8),
  C = pmin(100, hcl[,2] * 1.1),
  H = hcl[,3]))
# add some reds and yellows
cols_long <- c(cols_clean, plasma(20, begin = 0.3, end = 1), rocket(10, begin = 0.4, end = 0.8))
image(matrix(seq_along(cols_long), nrow = 1), col = cols_long, axes = FALSE)
# select final colors 
set.seed(1512) # map 1
#set.seed(15) # map 2
#set.seed(20) # map 3
cols_final <- sample(cols_long, n, replace=T)
image(matrix(seq_along(cols_final), nrow = 1), col = cols_final, axes = FALSE)

#_____________
# base map

finalPalette <- cols_final
names(finalPalette) <- unique(lines_sf$study_species) # assign names to know which color to use for inset
map <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey50", linewidth = 0.2) +
  geom_sf(data = lines_sf, aes(color = study_species), linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = finalPalette) +
  coord_sf(clip = "on") +
  theme_void() +
  theme(legend.position = "none")
map

# save for inkscape
library(svglite)
ggsave("Revision/NewSupplFigures/globalMap/newMap1_commuting.svg",
  plot = map, device = svglite, width = 10, height = 6)

#_____________
# Inset of selected species/dataset

unique(lines_inset$study_species)

# For turkey vulture
s <- "Cathartes aura MPIAB Cuba_Cathartes aura"
lines_s <- filter(lines_inset, study_species == s)
col <- finalPalette[names(finalPalette) == s]

# crop right extent
bb2 <- ext_to_bb(x=lines_s, f=0.2)
world_bb <- st_crop(world, bb2)

zoom2 <- ggplot() +
  geom_sf(data = world_bb, fill = "grey90", color = "grey50", linewidth = 0.2) +
  geom_sf(data = lines_s, color = col, linewidth = 0.5, alpha = 0.8) +
  geom_sf(data = bb2, fill = NA, color = "black", linewidth = 0.3, alpha = 0.8) +
  coord_sf(clip = "on") +
  theme_void() +
  theme(legend.position = "none")
zoom2

ggsave("Revision/NewSupplFigures/globalMap/insetTurkey.svg",
       plot = zoom2, device = svglite, width = 5, height = 2)


# For black legged kittiwake
s <- "AllisonPatterson_BLKI_Rissa tridactyla"
lines_s <- filter(lines_inset, study_species == s)
#commuting_s <- filter(commuting_inset, , study_species == s)
col <- finalPalette[names(finalPalette) == s]

# crop right extent
bb3 <- ext_to_bb(x=lines_s, f=0.2)
world_bb <- st_crop(world, bb3)

zoom3 <- ggplot() +
  geom_sf(data = world_bb, fill = "grey90", color = "grey50", linewidth = 0.2) +
  geom_sf(data = lines_s, color = col, linewidth = 0.5, alpha = 0.8) +
  geom_sf(data = bb3, fill = NA, color = "black", linewidth = 0.3, alpha = 0.8) +
  coord_sf(clip = "on") +
  theme_void() +
  theme(legend.position = "none")
zoom3

ggsave("Revision/NewSupplFigures/globalMap/insetKitty.svg",
       plot = zoom3, device = svglite, width = 5, height = 3)

# # For Sula
# s <- "AllisonPatterson_PEBO_Sula variegata"
# lines_s <- filter(lines_inset, study_species == s)
# segms_s <- filter(commuting_inset, study_species == s)
# 
# col <- finalPalette[names(finalPalette) == s]
# 
# # crop right extent
# bb1 <- ext_to_bb(x=lines_s, f=0.2)
# world_bb <- st_crop(world, bb1)
# 
# zoom_sula <- ggplot() +
#   geom_sf(data = world_bb, fill = "grey90", color = "grey50", linewidth = 0.2) +
#   geom_sf(data = lines_s, color = col, linewidth = 0.5, alpha = 0.8) +
#   #geom_sf(data = segms_s, color = col, linewidth = 0.5, alpha = 0.8) +
#   geom_sf(data = bb1, fill = NA, color = "black", linewidth = 0.3, alpha = 0.8) +
#   coord_sf(clip = "on") +
#   theme_void() +
#   theme(legend.position = "none")
# zoom_sula
# 
# ggsave("Revision/NewSupplFigures/globalMap/insetSula.svg",
#        plot = zoom_sula, device = svglite, width = 5, height = 6)



# Add box of the three insets to the world map
map_bbox <- map + 
  #geom_sf(data = bb1, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = bb2, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = bb3, fill = NA, color = "black", linewidth = 0.3)
map_bbox

ggsave("Revision/NewSupplFigures/globalMap/newMap1_commuting_insetsBB.svg",
       plot = map_bbox, device = svglite, width = 10, height = 6)

