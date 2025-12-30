

library(ggplot2)
library(dplyr)

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
#setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

speciesDf <- read.csv("DataFinalSummary/4supplementary_finalSpeciesList_wingMorphology_Feb2025_final.csv")
speciesDf$wingSpanAverage_birdsOfWorld_cm[speciesDf$species_phy=="Anser albifrons"] <- speciesDf$WingSpan_cm_MM2021[speciesDf$species_phy=="Anser albifrons"]
speciesDf$aspectRatio <- speciesDf$wingSpanAverage_birdsOfWorld_cm^2/speciesDf$wingArea_ellipse_cm2

# Add run swim categories from step 8 for coloring the plot
allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")
soarFlap <- group_by(allSegmDfs, species_phy) %>% summarise(soarFlap_pgls=unique(soarFlap_pgls))
speciesDf <- merge(speciesDf, soarFlap, by="species_phy", all.x=T)

# plot(speciesDf$aspectRatio, speciesDf$wingLoading_kgm2,
#      col = as.factor(speciesDf$soarFlap_pgls))
# text(speciesDf$aspectRatio, speciesDf$wingLoading_kgm2,
#      labels = speciesDf$species_phy, pos = 4, cex = 0.7)
# 
# plot(speciesDf$wingLoading_kgm2, speciesDf$Body_mass_kg,
#      col = as.factor(speciesDf$soarFlap_pgls))
# text(speciesDf$wingLoading_kgm2, speciesDf$Body_mass_kg,
#      labels = speciesDf$species_phy, pos = 4, cex = 0.7)
# 
# plot(speciesDf$wingArea_ellipse_cm2, speciesDf$Body_mass_kg,
#      col = as.factor(speciesDf$soarFlap_pgls))
# text(speciesDf$wingArea_ellipse_cm2, speciesDf$Body_mass_kg,
#      labels = speciesDf$species_phy, pos = 4, cex = 0.7)
# 
# plot(speciesDf$aspectRatio, speciesDf$Body_mass_kg,
#      col = as.factor(speciesDf$soarFlap_pgls))
# text(speciesDf$aspectRatio, speciesDf$Body_mass_kg,
#      labels = speciesDf$species_phy, pos = 4, cex = 0.7)

# Plotting the same but on a log curve
hist(speciesDf$avg_probFlap)
hist(speciesDf$Body_mass_kg)
m <- lm(log(avg_probFlap)~log(Body_mass_kg), data=speciesDf)
summary(m)
plot(m)
coef(m)

# now model only the lowest 1st percentile to draw this lower limit 
# of the lowest possible prop flapping per body mass
# Fit the quantile regression model for the 1st percentile
library(quantreg)
modLowerLim <- rq(log(avg_probFlap)~log(Body_mass_kg), data=speciesDf, tau = 0.01)
summary(modLowerLim)
text(0, -2, labels = paste("y=", round(coef(modLowerLim)[2], 2), "x", round(coef(modLowerLim)[1], 2), sep=""), pos = 4, cex = 1.5, col = "firebrick")


# New figure adding runSwim colors and adding upper and lower limit in proportion of flapping
library(Cairo)
myCols <- c("gold2","dodgerblue3")

speciesDf$species_abb <- paste0(substr(speciesDf$species,1,1), ". ",
                                   sapply(strsplit(as.character(speciesDf$species), " "), "[",2))

#CairoPDF("Plots/finalPlots/newModelPlot_newRunSwim_phyloModels/bodyMassVSpropFlap_loglog_newLowerLimitFlap.pdf", width=6, height=6)
CairoPDF("Plots/finalPlots/newModelPlot_newRunSwim_phyloModels/bodyMassVSpropFlap_loglog_newLowerLimitFlap_noSpeciesNames.pdf", width=6, height=6)
ggplot(speciesDf, mapping=aes(x=log(Body_mass_kg), y=log(avg_probFlap), color=runSwim, label=species_abb)) +
  xlim(-1.6,3) + ylim(-6,0.5) +
  geom_abline(slope = 0, intercept = 0, col="gold", linetype = 2, size = 1, alpha=0.7) +
  geom_abline(slope = coef(modLowerLim)[2], intercept = coef(modLowerLim)[1], col="dodgerblue1", linetype = 2, size = 1, alpha=0.7) +
  geom_point(alpha=1, size=3.5, col="white") +
  geom_point(alpha=1, size=2.5) +
  #geom_text(aes(label=species_abb), hjust=-0, vjust=-1, col="black", size=3.5, fontface = "italic") +
  theme_bw() + xlab("\nBody mass (kg)") + ylab("Probability of flapping flight\n") +
  scale_color_manual(values=myCols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none")
dev.off()


#___________________________________________________________________________
## Plot the distribution of all trajectories of all downloaded studies ####
#___________________________________________________________________________

setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")

library(ggplot2)
library(lubridate)
library(data.table)
theme_set(theme_bw())
library(pals) # palette library for class data
library(colorspace)
library(sf)
library(maps)
world_map <- map_data("world")

# Import final dataset
finalDF <- readRDS("DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV_Feb2025.rds")
sub <- finalDF[,c(1:15)]

# Create a pastel palette
rainb_col <- rainbow(length(unique(sub$species)))
rainb_pastel <- lighten(desaturate(rainb_col, 0.3), 0.1)
pie(rep(1, 38), col = rainb_pastel, labels = 1:38) #visualise it

myCols <- as.vector(rainb_pastel)

#pdf("Plots/finalPlots/distrMaps_allStudies_gps.pdf", 10, 6.5)
tiff("Plots/finalPlots/distrMaps_allStudies_gps_grey.tif", width=10, height=6.5, res=400, unit="in")
ggplot() +
  #geom_polygon(world_map, mapping = aes(x = long, y = lat, group = group), fill="grey10", colour = "darkgrey") +
  geom_polygon(world_map, mapping = aes(x = long, y = lat, group = group), fill="grey70", colour = "grey35") +
  #geom_polygon(world_map, mapping = aes(x = long, y = lat, group = group), fill="beige", colour = "grey70") +
  geom_point(data=sub, mapping = aes(x=location.long, y=location.lat, color=species, group=study.id), alpha=0.3, size=0.4) +
  theme_bw() + xlab("") + ylab("") + ylim(c(-65,90)) +
  scale_color_manual(values=myCols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #legend.keys.align=0.5,
        legend.position="bottom", legend.direction="vertical", legend.title.align=0.5,
        legend.title = element_text(size=8), legend.text=element_text(size=6)
        #,panel.background = element_rect(fill = "beige")
  ) +
  guides(col=guide_legend("Species", ncol=8, override.aes = list(size=2))) #change point size only in legend
dev.off()

#______________________________________________________________
# Zoom in to one individual of one species for composite figure

oneInd <- sub[sub$study.name=="Demoiselle Crane High Resolution Mongolia" & 
                sub$individual.local.identifier=="crane 6576" & month(sub$timestamp)%in%c(2:3),] #1:3 or 10:12
bbox <- c(range(oneInd$location.long),range(oneInd$location.lat))
colIdx <- which(levels(as.factor(sub$species)) == unique(oneInd$species))
asp_ratio <- diff(bbox[3:4]) / diff(bbox[1:2])  # height / width

ggplot(data=oneInd, aes(x=location.long, y=location.lat), alpha=0) +
  geom_polygon(world_map, mapping = aes(x = long, y = lat, group = group), fill="grey70", colour = "grey35") +
  geom_path(data=oneInd, mapping = aes(x=location.long, y=location.lat, group=study.id), col=myCols[colIdx], size=0.4) +
  geom_point(data=oneInd, mapping = aes(x=location.long, y=location.lat, group=study.id), col=myCols[colIdx], fill=myCols[colIdx], shape=21, size=0.4) +
  theme_bw() + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #legend.keys.align=0.5,
        legend.position="none") +
  coord_fixed(xlim = bbox[1:2], ylim = bbox[3:4])
ggsave("Plots/finalPlots/distrMaps_allStudies_gps_grey_zoomInCrane.pdf", width = 2, height = 2 * asp_ratio, units = "in")

#__________________________
## Explore the Osprey ####
#__________________________

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")

geo_mean <- function(x) exp(mean(log(x), na.rm = TRUE))

osp <- finalDF[finalDF$species=="Pandion haliaetus",-c(40,42)]
summary(osp$height.above.msl)

meanPerSegm <- osp %>% group_by(track_flight_id, individual.local.identifier) %>% 
  summarise(flapping_prob=mean(flapping_prob), groundSpeed_ms=mean(groundSpeed_ms), nPoints=n()) %>% ungroup()
meanPerInd <- meanPerSegm %>% group_by(individual.local.identifier) %>% 
  summarise(flapping_prob=mean(flapping_prob), groundSpeed_ms=mean(groundSpeed_ms), nPoints=sum(nPoints))
mean(meanPerInd$flapping_prob)
geo_mean(meanPerInd$flapping_prob)
mean(log(meanPerInd$flapping_prob))

bbox <- c(range(osp$location.long),range(osp$location.lat))
asp_ratio <- diff(bbox[3:4]) / diff(bbox[1:2])  # height / width
ospCols <- rainbow(length(unique(osp$individual.local.identifier)))
ospCols <- lighten(desaturate(ospCols, 0.3), 0.1)

ggplot() +
  #geom_polygon(world_map, mapping = aes(x = long, y = lat, group = group), fill="grey70", colour = "grey35") + #as dataframe
  geom_sf(data = world, fill = "grey70", color = "grey35") + # as sf
  geom_point(data=osp, mapping = aes(x=location.long, y=location.lat, fill=individual.local.identifier, col=individual.local.identifier), shape=21, size=0.4) +
  scale_color_manual(values=ospCols) +
  theme_bw() + xlab("") + ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #legend.keys.align=0.5,
        legend.position="none") +
  coord_sf(xlim = bbox[1:2], ylim = bbox[3:4], expand=F)

# Transform points in sf 
osp_sf <- st_as_sf(osp, coords = c("location.long", "location.lat"), crs=st_crs(world))

# calculate proportion of points per individual over sea vs over land
world_union <- st_union(world, by_feature=F)
ggplot() + geom_sf(data=world_union, fill = "grey70")

# Calculate the proportion of points outside the polygon (sea) versus land
# Per segment
osp_sf_segm <- split(osp_sf, osp_sf$track_flight_id)
atSea_perSegm <- rbindlist(lapply(osp_sf_segm, function(x){
  inside <- st_within(x, world) #elements == 0 are outside
  n_outside <- sum(lengths(inside) == 0)
  n_total <- nrow(x)
  prop_outside <- n_outside / n_total # Proportion outside
  return(data.frame(track_flight_id = unique(x$track_flight_id),
                    prop_overSea = prop_outside))
}))

meanPerSegm_sea <- merge(meanPerSegm, atSea_perSegm, by="track_flight_id")
meanPerSegm_sea
plot(flapping_prob~prop_overSea, meanPerSegm_sea[-which(meanPerSegm_sea$individual.local.identifier=="IBH_Infiernillo-181171_181171"),])
plot(groundSpeed_ms~prop_overSea, meanPerSegm_sea[-which(meanPerSegm_sea$individual.local.identifier=="IBH_Infiernillo-181171_181171"),])

# Per individual
osp_sf_ind <- split(osp_sf, osp_sf$individual.local.identifier)
atSea_perInd <- rbindlist(lapply(osp_sf_ind, function(x){
  inside <- st_within(x, world) #elements == 0 are outside
  n_outside <- sum(lengths(inside) == 0)
  n_total <- nrow(x)
  prop_outside <- n_outside / n_total # Proportion outside
  return(data.frame(individual.local.identifier = unique(x$individual.local.identifier),
                    prop_overSea = prop_outside))
}))


meanPerInd_sea <- merge(meanPerInd, atSea_perInd, by="individual.local.identifier")
meanPerInd_sea
plot(flapping_prob~prop_overSea, meanPerInd_sea[-8,])
plot(groundSpeed_ms~prop_overSea, meanPerInd_sea[-8,])

# test averages per segment
geo_mean <- function(x) exp(mean(log(x), na.rm = TRUE))
osp <- allSegmDfs[allSegmDfs$species_phy=="Pandion haliaetus",]
perInd <- osp %>% group_by(individualID) %>% 
  summarise(gmean_pFlap=geo_mean(avg_probFlap),mean_pFlap=mean(avg_probFlap, na.rm=T),
            totNlocs=sum(totNlocs), duration=sum(segm_timeDuration_h)) %>% 
  as.data.frame()
geo_mean(perInd$gmean_pFlap);mean(perInd$mean_pFlap, na.rm=T) # per individual and species
geo_mean(perInd$gmean_pFlap[-8]);mean(perInd$mean_pFlap[-8], na.rm=T) # per individual and species
geo_mean(osp$avg_probFlap);mean(osp$avg_probFlap, na.rm=T) # per species

perInd$gmean_pFlap <- round(perInd$gmean_pFlap,3)
perInd$mean_pFlap <- round(perInd$mean_pFlap,3)
perInd

#_____________________________________
# Summary values/tables for manuscript
#_____________________________________

library(dplyr)
geo_mean <- function(x) exp(mean(log(x), na.rm = TRUE))

# Import data without outlier individuals
allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025_noOutliers.rds")

summary(allSegmDfs$obs_tot_COT_J_KgM)

#___________________
# Median residuals

allSegmDfs %>% group_by(soarFlap_pgls) %>% summarise(mean(COT_residuals), median(COT_residuals))

#______________________________________
# Summary per individual and species
allSegmDfs$log_probFlap <- log(allSegmDfs$avg_probFlap)
MRperID <- allSegmDfs %>% group_by(individualID) %>%
  summarise(species_phy = unique(species_phy),
            obs_MR_pass_W=mean(obs_MR_pass_W, na.rm=T),
            obs_MR_flap_W=mean(obs_MR_flap_W, na.rm=T),
            tot_MR_W=mean(tot_MR_W, na.rm=T),
            obs_tot_COT_J_KgM=mean(obs_tot_COT_J_KgM, na.rm=T),
            geomMean_COT=mean(log(obs_tot_COT_J_KgM), na.rm=T),
            avg_grSpeed_ms=mean(avg_grSpeed_ms, na.rm=T),
            avg_probFlap=mean(avg_probFlap, na.rm=T),
            geomMean_probFlap=mean(log_probFlap, na.rm=T),
            soarFlap_pgls=unique(soarFlap_pgls, na.rm=T),
            Body_mass_kg=unique(Body_mass_kg))
nrow(MRperID)

MRperSpecies <- MRperID %>% group_by(species_phy) %>%
  summarise(obs_MR_pass_W=mean(obs_MR_pass_W, na.rm=T),
            obs_MR_flap_W=mean(obs_MR_flap_W, na.rm=T),
            tot_MR_W=mean(tot_MR_W, na.rm=T),
            obs_tot_COT_J_KgM=mean(obs_tot_COT_J_KgM, na.rm=T),
            avg_grSpeed_ms=mean(avg_grSpeed_ms, na.rm=T),
            avg_probFlap=mean(avg_probFlap, na.rm=T),
            geomMean_probFlap=mean(geomMean_probFlap, na.rm=T),
            geomMean_COT=mean(geomMean_COT, na.rm=T),
            soarFlap_pgls=unique(soarFlap_pgls, na.rm=T),
            Body_mass_kg=unique(Body_mass_kg))

MRperSpecies %>% arrange(obs_tot_COT_J_KgM) %>% print(n=Inf, width=Inf)
MRperSpecies$taxon <- "Birds"
MRperSpecies$taxon[grepl("Eidolon|Pteropus", MRperSpecies$species_phy)] <- "Bats"
MRperSpecies <- MRperSpecies %>% arrange(desc(taxon), as.character(species_phy))

write.csv(MRperSpecies, file="DataFinalSummary/COT_summaryperSpecies.csv", row.names = F)

#__________________________________________________________
# Summary per species separately for flappers and soarers

COTperGroup <- allSegmDfs %>% group_by(species_phy, soarFlap_pgls) %>%
  summarise(avg_obsCOT=mean(obs_tot_COT_J_KgM),
            IQR_obsCOT=IQR(obs_tot_COT_J_KgM),
            min_obsCOT=min(obs_tot_COT_J_KgM),
            max_obsCOT=max(obs_tot_COT_J_KgM),
            avg_probFlap=mean(avg_probFlap),
            min_probFlap=min(avg_probFlap),
            max_probFlap=max(avg_probFlap),
            mean_COTresid=mean(COT_residuals),
            wingLoad=unique(wingLoading_kgm2),
            nSegments=n(),
            bodyMass_kg=unique(Body_mass_kg))

# check highest and lowest prop flapping
COTperGroup %>% arrange(avg_probFlap) %>% select(species_phy,avg_probFlap) %>% print(n=Inf)

# check highest and lowest COT residuals
COTperGroup %>% arrange(mean_COTresid) %>% select(species_phy,mean_COTresid) %>% print(n=Inf, width=Inf)

# check mean wing loading per group
COTperGroup %>% group_by(soarFlap_pgls) %>% summarise(mean(wingLoad, na.rm=T))

# check highest and lowest COT
COTperGroup %>% arrange(avg_obsCOT) %>% print(n=Inf, width=Inf)

# check IQR of obsCOT
COTperGroup[which.max(COTperGroup$IQR_obsCOT),]
COTperGroup[which.min(COTperGroup$IQR_obsCOT),]
COTperGroup %>% group_by(soarFlap_pgls) %>% summarise(mean(IQR_obsCOT))
COTperGroup %>% arrange(soarFlap_pgls,IQR_obsCOT) %>% select(soarFlap_pgls,IQR_obsCOT,wingLoad) %>% print(n=Inf)

# Relationship between IQR of COT and body mass
plot(IQR_obsCOT~bodyMass_kg, data=COTperGroup)
plot(log(IQR_obsCOT)~log(bodyMass_kg), COTperGroup, 
     col=c("gold2","dodgerblue2")[as.factor(COTperGroup$soarFlap_pgls)], pch=20)

modIQR1 <- lm(log(IQR_obsCOT)~log(bodyMass_kg)*soarFlap_pgls,COTperGroup)
summary(modIQR1)

# Relationship between IQR of COT and wing loading
plot(IQR_obsCOT~wingLoad, data=COTperGroup)
plot(log(IQR_obsCOT)~log(wingLoad),COTperGroup,
     col=c("gold2","dodgerblue2")[as.factor(COTperGroup$soarFlap_pgls)], pch=20)

modIQR2 <- lm(log(IQR_obsCOT)~log(wingLoad)*soarFlap_pgls,COTperGroup)
summary(modIQR2)


#_______________________________________________________________
# Compare swan and condor, speed and daily distance travelled:
cyg <- allSegmDfs[allSegmDfs$species_phy=="Cygnus buccinator",]; cyg <- cyg[date(cyg$segm_timeStart) == date(cyg$segm_timeEnd),]
con <- allSegmDfs[allSegmDfs$species_phy=="Vultur gryphus",]
summary(cyg[,c("segm_timeDuration_h","avg_stepLength_m","avg_grSpeed_ms","tot_distCovered_km","StDev_grSpeed_ms","wingLoading_kgm2","Body_mass_kg")])
summary(con[,c("segm_timeDuration_h","avg_stepLength_m","avg_grSpeed_ms","tot_distCovered_km","StDev_grSpeed_ms","wingLoading_kgm2","Body_mass_kg")])

table(date(cyg$segm_timeStart) == date(cyg$segm_timeEnd))
table(date(con$segm_timeStart) == date(con$segm_timeEnd))

cygSum <- cyg %>% group_by(individualID, date(segm_timeStart)) %>%
  summarise(dailyDist_km=sum(tot_distCovered_km),
            dailySpeed=mean(avg_grSpeed_ms)) %>% ungroup() %>%
  summarise(minDist=min(dailyDist_km), maxDist=max(dailyDist_km),medDist=median(dailyDist_km),
            minSpeed=min(dailySpeed), maxSpeed=max(dailySpeed),medSpeed=median(dailySpeed))
conSum <- con %>% group_by(individualID, date(segm_timeStart)) %>%
  summarise(dailyDist_km=sum(tot_distCovered_km),
            dailySpeed=mean(avg_grSpeed_ms)) %>% ungroup() %>%
  summarise(minDist=min(dailyDist_km), maxDist=max(dailyDist_km),medDist=median(dailyDist_km),
            minSpeed=min(dailySpeed), maxSpeed=max(dailySpeed),medSpeed=median(dailySpeed))
cygSum
conSum

