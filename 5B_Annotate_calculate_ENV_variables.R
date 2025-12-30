
#___________________
# Annotate data ####
#___________________

setwd("/media/mscacco/Cyprus")
#cot_path <- "/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT/" #laptop
cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/" #desktop
source(paste0(cot_path,"Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R")) # functions also include necessary libraries

dir.create(paste0(cot_path,"DataFinalSummary/AnnotatedData"))

# Load the file to annotate
finalDf <- readRDS(paste0(cot_path,"DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs.rds"))
length(unique(finalDf$mergingCol)) == nrow(finalDf) # mergingCol is the unique row ID that we will use for the extraction
length(unique(finalDf$event.id)) == nrow(finalDf)
finalDf$event.id <- NULL

#_____________________
# Annotate topography

topoMerge <- rast(paste0(cot_path,"DEM/world_DEM_slope_aspect_deFerranti_1km.tif")) 
dem <- extract(topoMerge[["dem"]], finalDf[,c("location.long","location.lat")], method="bilinear")
aspSlope <- extract(topoMerge[[c("slope","aspect")]], finalDf[,c("location.long","location.lat")], method="simple")
nrow(finalDf)==nrow(aspSlope)
finalDf <- cbind(finalDf, dem, aspSlope)

summary(finalDf$dem)
summary(finalDf$aspect)

saveRDS(finalDf, paste0(cot_path,"DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV.rds"))

# Check which species (and where) have NAs in the topography (it is animals flying over water)
table(finalDf$species[is.na(finalDf$slope)])
plot(topoMerge[[1]])
points(finalDf$location.long[is.na(finalDf$slope)], finalDf$location.lat[is.na(finalDf$slope)], col="red", pch=19, cex=0.5)
plot(topoMerge[[1]], xlim=c(0,100), ylim=c(50,90))
points(finalDf$location.long[is.na(finalDf$slope)], finalDf$location.lat[is.na(finalDf$slope)], col="red", pch=19, cex=0.5)


#__________________________
# Define height to extract

finalDf <- readRDS(paste0(cot_path,"DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV.rds"))

# Define height to download
# Let's extract the median height above ellipsoid/sea level just to decide on a pressure level to download wind data from:
summary(finalDf$height_gener) #median == 1245 m ~= 872.68 mbar let's use 1000 m ~= 898.75 mbar
# For the conversion we used: https://www.herramientasingenieria.com/onlinecalc/altitude/altitude.html
#Pressure levels in ECMWF ERA-5: 
#1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50/30/20/10/7/5/3/2/1
finalDf$annotationHeightEll_m <- 1000

# Set the directory path for the ERA5 files
# List ERA5 files in one of the directory, (files dates are the same for single levels and pressure levels)
era5_files <- list.files("ERA5downloads_singleLev", pattern = "^download_ERA5_.*singleLev.nc$", full.names = TRUE)

#__________________________________________________________________
# Extract track data for annotation, for each ERA5 file in parallel

registerDoMC(cores = parallel::detectCores() - 3) # Register a parallel back-end

start_time <- Sys.time()
extracted_data_list <- llply(era5_files, function(x) {
  extract_track_and_era5_data(track=finalDf, 
                              era5_file=x,
                              eventCol = "mergingCol",
                              timeCol = "timestamp", 
                              coordsCol = c("location.long","location.lat"), 
                              heightCol = "annotationHeightEll_m")
}, .parallel = TRUE)
end_time <- Sys.time()
print(paste("Time taken for extracting track data:", round(difftime(end_time, start_time, units="mins"), 2), "minutes"))

# For every timestamp, we downloaded era5 data for the hour before and after to make time interpolation possible
# This means that we have some era5 data that do not have a corresponding section in the trajectory?
# Or why do we have elements with 0 rows in track data?
# Remove NULL elements from the list
extracted_data_list <- extracted_data_list[!sapply(extracted_data_list, is.null)]
# Also do not keep list elements with empty data frames
extracted_data_list <- extracted_data_list[sapply(extracted_data_list, function(x) nrow(x$track_data) > 0)]
# We make sure all rows in the dataset are associated to the era5 files
sum(sapply(lapply(extracted_single, "[[", 4), nrow)) == nrow(finalDf)

# We only need to do the extraction once, and just replace the file name, as the dates for singleLev and pressLev are exactly the same for both list of files
s <- list.files("ERA5downloads_singleLev", pattern = "^download_ERA5_.*singleLev.nc$", full.names = TRUE)
p <- list.files("ERA5downloads_pressLev", pattern = "^download_ERA5_.*pressLev.nc$", full.names = TRUE)
w <- list.files("ERA5downloads_upliftProxies", pattern = "^download_ERA5_.*star.nc$", full.names = TRUE)

sD <- sapply(strsplit(s, "_"), "[", 4)
pD <- sapply(strsplit(p, "_"), "[", 4)
wD <- sapply(strsplit(w, "_"), "[", 4)
identical(sD, pD)
identical(sD, wD)

# Duplicate the extraction list, creating one for the press lev and one for the single lev by changing the file names
extracted_single <- extracted_data_list

extracted_uplifts <- lapply(extracted_data_list, function(x){
  x$era5_file_name <- sub("singleLev","upliftProxies",x$era5_file_name)
  x$era5_file_name <- sub("singleLev","Wstar",x$era5_file_name)
  x$dataset <- "singleLev" # uplift proxies work as single level variables
  return(x)
})

extracted_press <- lapply(extracted_data_list, function(x){
  x$era5_file_name <- gsub("singleLev","pressLev",x$era5_file_name)
  x$dataset <- "pressLev" # pressure level variables work different as they get vertically interpolated
  return(x)
})

save(extracted_single, extracted_uplifts, extracted_press, file = paste0(cot_path,"DataFinalSummary/AnnotatedData/list_extractedFiles_to_annotate.rdata"))


#______________________________________
# Annotate ECMWF variables in parallel

setwd("/media/mscacco/Cyprus")
cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/" #desktop
source(paste0(cot_path,"Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R")) # functions also include necessary libraries

load(paste0(cot_path,"DataFinalSummary/AnnotatedData/list_extractedFiles_to_annotate.rdata")) #extracted_single, extracted_uplifts, extracted_press

registerDoMC(cores = parallel::detectCores() - 2) # Register a parallel back-end

## Single levels:
start_time <- Sys.time()
era5_data_single <- rbindlist(llply(extracted_single, 
                                    annotate_era5_data, 
                                    #pathToFolder = "DataFinalSummary/AnnotatedData", #only for really big datasets
                                    .parallel = TRUE))
# to check for errors run
#test <- lapply(extracted_single[1:10],function(x)try({annotate_era5_data(x)}))
end_time <- Sys.time()
print(paste("Time taken for annotating era5 data:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
summary(era5_data_single)
head(era5_data_single)
saveRDS(era5_data_single, file = paste0(cot_path,"DataFinalSummary/AnnotatedData/annotatedDf_era5_singleLev.rds"))

## Uplift proxy:
start_time <- Sys.time()
era5_data_uplifts <- rbindlist(llply(extracted_uplifts, 
                                     annotate_era5_data, 
                                     .parallel = TRUE))
end_time <- Sys.time()
print(paste("Time taken for annotating era5 data:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
summary(era5_data_uplifts)
head(era5_data_uplifts)
saveRDS(era5_data_uplifts, file = paste0(cot_path,"DataFinalSummary/AnnotatedData/annotatedDf_era5_Wstar.rds")) 

## Pressure levels:
start_time <- Sys.time()
era5_data_press <- rbindlist(llply(extracted_press, 
                                   annotate_era5_data, 
                                   .parallel = TRUE))
end_time <- Sys.time()
print(paste("Time taken for annotating era5 data:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
summary(era5_data_press)
head(era5_data_press)
saveRDS(era5_data_press, file = paste0(cot_path,"DataFinalSummary/AnnotatedData/annotatedDf_era5_pressLev.rds"))

#______________________________________
# Merge track info with era5 data ####
#______________________________________

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT/")
source("Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R") # functions also include necessary libraries

# Load the file to annotate
finalDf <- readRDS("DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV.rds")

# Load annotated files
singleLev <- readRDS("DataFinalSummary/AnnotatedData/annotatedDf_era5_singleLev.rds")
pressLev <- readRDS("DataFinalSummary/AnnotatedData/annotatedDf_era5_pressLev.rds")
Wstar <- readRDS("DataFinalSummary/AnnotatedData/annotatedDf_era5_Wstar.rds")
nrow(singleLev); nrow(pressLev); nrow(Wstar)

# Merge with original dataset
finalDf_era5 <- merge(finalDf, singleLev[,-c(2:5)], by.x="mergingCol", by.y="event_id", all.x=T)
finalDf_era5 <- merge(finalDf_era5, pressLev[,-c(2:5)], by.x="mergingCol", by.y="event_id", all.x=T)
finalDf_era5 <- merge(finalDf_era5, Wstar[,-c(2:5)], by.x="mergingCol", by.y="event_id", all.x=T)
head(finalDf_era5)

# Calculate orographic uplift potential from topo and wind columns
finalDf_era5$Orographic_uplift_potential <- oroUplift(u = finalDf_era5$u10_10_metre_U_wind_component,
                                                      v = finalDf_era5$v10_10_metre_V_wind_component,
                                                      slope = (finalDf_era5$slope * pi / 180), # in rad
                                                      aspect = (finalDf_era5$aspect * pi / 180)) # in rad between 0 and 2pi
hist(finalDf_era5$Orographic_uplift_potential)
summary(finalDf_era5)

# Add also wind speed and wind direction
finalDf_era5$windSpeed_1000m <- windSpeed(u = finalDf_era5$u_U_component_of_wind,
                                          v = finalDf_era5$v_V_component_of_wind)
finalDf_era5$windDirFROM_1000m <- windDir_from_deg(u = finalDf_era5$u_U_component_of_wind,
                                                   v = finalDf_era5$v_V_component_of_wind)

# Calculate wind support, crosswind and airspeed
finalDf_era5$windSupport_1000m <- windSupport(u = finalDf_era5$u_U_component_of_wind,
                                              v = finalDf_era5$v_V_component_of_wind,
                                              dg = finalDf_era5$segmentDir)

finalDf_era5$crossWind_1000m <- crossWind(u = finalDf_era5$u_U_component_of_wind,
                                          v = finalDf_era5$v_V_component_of_wind,
                                          dg = finalDf_era5$segmentDir)

finalDf_era5$airspeed_1000m <- airspeed(Vg = finalDf_era5$groundSpeed_ms, 
                                        Ws = finalDf_era5$windSupport_1000m, 
                                        Cw = finalDf_era5$crossWind_1000m)

# Some checks
summary(finalDf_era5[,c(41,43,44,57:62)])
boxplot(Orographic_uplift_potential~flapping_bin, data=finalDf_era5)
boxplot(`w_star_Thermal_uplift_potential_(W*)`~flapping_bin, data=finalDf_era5)

# Save final dataset for the models in next steps
#saveRDS(finalDf_era5, "DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV.rds")
#saveRDS(finalDf_era5, "DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV_Dec2024.rds")
saveRDS(finalDf_era5, "DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs_ENV_Feb2025.rds")

