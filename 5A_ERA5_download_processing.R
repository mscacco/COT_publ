## Step 4B. Download and annotate ERA5 data
## Script and functions modified from: https://github.com/kamransafi/ERA5_move2
## We will download and annotate U_wind, V_wind, and all variables needed for the calculation of oroUplift and thermUplift
## To run this script the dataset resulting from the script "3B.VEDBAclassification_MixR_createFinalDf.R" is needed

#________________________
# Download ERA5 data ####
#________________________

setwd("/media/mscacco/Cyprus")
#cot_path <- "/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT" #laptop
cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/" #desktop
source(paste0(cot_path,"Scripts/COT/Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R"))

finalDf <- readRDS(paste0(cot_path,"DataFinalSummary/FinalDf_March2024_perPoint_VedbaGs_filteredSpecies_flappingProbs.rds"))
dir.create("ERA5downloads")
localPath <- "ERA5downloads"

# Store credentials, you can find them (or request them) here: https://cds.climate.copernicus.eu/user/login?destination=%2F%23!%2Fhome
wf_set_key(user = "",
           key = "",
           service = "cds")



# Create a request table, specifying the dates and hours and the area of each request
# In this case we want monthly request for the entire area (bbox) covered by our dataset
# the function takes either a move2 object or a data.frame (for a data.frame, the column names of time and coordinates need to be specified unless they fit the default of a Movebank dataframe)
reqTable <- timeSpace2request(finalDf, timeUnit = "day", area = "bbox",
                              timeCol = "timestamp",
                              coordsCol = c("location.long","location.lat"))
print(paste("Number of days to be requested:", nrow(reqTable)))

# The var geopotential should always be downloaded among the press lev variables for correct estimation of the pressure level height
vars_pressLev <- c("Geopotential", "Specific humidity", "Temperature", "u_component_of_wind" , "v_component_of_wind")
levels <- c("650", "750", "850", "900", "950", "1000")
vars_singleLev <- c("2m_temperature", "boundary_layer_height", "instantaneous_moisture_flux",
                    "instantaneous_surface_sensible_heat_flux", "surface_pressure",
                    "10m_u_component_of_wind", "10m_v_component_of_wind")

# Create two separate requests for the surface variables and the pressure levels variables
# Each element of each list contains monthly requests
# reqList_pressLev <- create_request_list(reqTable, 
#                                         datasetName = "reanalysis-era5-pressure-levels",
#                                         vars = vars_pressLev, 
#                                         levels = levels)
reqList_singleLev <- create_request_list(reqTable, 
                                         datasetName = "reanalysis-era5-single-levels",
                                         vars = vars_singleLev)

# remove partial downloads
unlink(list.files(localPath, pattern="ecmwf", full.names = TRUE))

# Check for existing files and remove corresponding existing requests from the request list
# For press levels:
# existing_files <- list.files(localPath, pattern = "^download_ERA5_.*\\pressLev.nc$", full.names = TRUE)
# reqList_pressLev <- reqList_pressLev[!vapply(reqList_pressLev, function(req) file.path(localPath, req$target) %in% existing_files, FUN.VALUE = logical(1))]
# print(paste("Number of total requests to be submitted:", length(reqList_pressLev)))
# For single levels:
existing_files <- list.files(localPath, pattern = "singleLev.nc", full.names = TRUE)
reqList_singleLev <- reqList_singleLev[!vapply(reqList_singleLev, function(req) file.path(localPath, req$target) %in% existing_files, FUN.VALUE = logical(1))]
print(paste("Number of total requests to be submitted:", length(reqList_singleLev)))

# Adjust the timeout option of curl in case of internet connection problems
# options(timeout = 20000)
# library(httr)
# set_config(config(timeout = 300))
# library(curl)
# handle <- new_handle() ## curl::curl_options()
# handle_setopt(handle, .list = list(connecttimeout = 300, timeout = 900))

# Submit requests to CDS API based on the list of requests
# send_requests(reqList_pressLev, localPath, user=CDS_USER_ID,
#               batch_size = 5, requestFileName = "all_requests_pressLev")
send_requests(reqList_singleLev, localPath, user=CDS_USER_ID,
              batch_size = 5, requestFileName = "all_requests_singLev")

# Login here to check the status of the requests:
# https://cds.climate.copernicus.eu/cdsapp#!/yourrequests


#_____________________________________
# Calculate W* for each timestep ####
#____________________________________

#_________________________________________________________________________
# Add W* calculation as additional variable in the single level variables

setwd("/media/mscacco/Cyprus")
cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/" #desktop
source(paste0(cot_path,"Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R"))

# This code assumes that the surface level variables and the pressure level variables were downloaded for the exact same dates and time steps
# It also assumes that the following variables were downloaded for the two groups:
# vname_surface <- c("t2m", "sp", "ishf", "ie", "blh") #"2 metre temperature","Surface pressure","Instantaneous surface sensible heat flux","Instantaneous moisture flux","Boundary layer height" 
# vname_pressure <- c("t", "q", "z") # "Temperature", "Specific humidity", "Geopotential"
# Note, the pressure levels are in hPa == mbar
# Surface pressure is in Pascals, so Pa/100 == mbar
# Geopotential height in m == Geopotential/gravitational acceleration

# For the same dates we need both surface level and pressure level variables, so we extract them by date
uniqueDates <- sapply(strsplit(list.files("ERA5downloads_singleLev", pattern="nc"), "_"), "[", 3)

# Create a new folder to store the w_star. 
dir.create("ERA5downloads_upliftProxies")
upliftPath <- "ERA5downloads_upliftProxies"
# For simplicity W* will be added as an extra variable in the rasters including the single layer variables
# dir.create("ERA5downloads_singleLev_Wstar") # do not change this folder name (it is used for the output file in the lapply below)

surfFls <- list.files("ERA5downloads_singleLev", pattern = "nc", full.names = T)
pressFls <- list.files("ERA5downloads_pressLev", pattern = "nc", full.names = T)

# For each day, we create a raster with hourly W* per pixel (as many layers as the number of hours)
library(doParallel)
library(plyr)
registerDoParallel(6)

llply(uniqueDates, function(date){ # per date
  
  # Load rasters for each date
  surf <- rast(grep(date, surfFls, value=T))
  press <- rast(grep(date, pressFls, value=T))
  # Split by hour/time-step (both press and surf variables were downloaded for the same hours/time-steps)
  surf_ls <- split(surf, time(surf))
  press_ls <- split(press, time(press))
  
  # Each level has a constant pressure, but the corresponding height in m vary depending on the Geopotential (gravitational potential energy of a unit mass at a particular location m2/s2, relative to mean sea level). 
  # To obtain geopotential height in m, divide the Geopotential (in m2/s2) by gravitational acceleration (9.80665 m/s2) https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
  # From this lapply we will obtain a new surface raster stack for the above date, with one element per timestep, and an added layer corresponding to w_star per tiemstep 
  wstar_ls <- lapply(seq_along(surf_ls), function(i){ # per hour/time step
    
    if(unique(time(surf_ls[[i]])) != unique(time(press_ls[[i]]))){stop("The time step of surface variable and pressure variable does not match!")}
    
    geop_r <- press_ls[[i]][[grep("z_", names(press_ls[[i]]))]] # extract raster of geopotential, in m2/s2
    
    pressLevel_height_m <- geop_r / 9.80665 # Geopotential height, in m, which is the height in m corresponding to each pressure level
    # extract values of blh and pressure level vars. The values are returned as a vector/matrix, in cell-order by layer. 
    blh <- values(surf_ls[[i]][[grep("blh", names(surf_ls[[i]]))]]) # Boundary layer height, in m
    # these are per pressure level
    geop_df <- values(geop_r) # Geopotential, in m2/s2
    hum_df <- values(press_ls[[i]][[grep("q_", names(press_ls[[i]]))]]) # Specific humidity, in Kg/Kg (expressed as Kg of water vapor present in each kilogram of air)
    temp_df <- values(press_ls[[i]][[grep("t_", names(press_ls[[i]]))]]) # Temperature, in K
    # exctract vars at pressure level height that is closest to the blh in m
    pressLevel_height_df <- values(pressLevel_height_m) # make a matrix of heights per pressure level
    colnames(pressLevel_height_df) <- gsub("_.*","",gsub(".*=","",colnames(pressLevel_height_df))) # name column with pressure level value
    # take the absolute difference of each height value from the bhl
    diff_df <- abs(sweep(pressLevel_height_df, 1, blh, FUN="-"))
    # take the minimum difference per row and assign the corresponding pressure level
    whichCol <- apply(diff_df, 1, which.min)
    
    # extract the press level closest to bhl in hPa
    pressLev <- as.numeric(colnames(diff_df)[whichCol]) # from hPa == mbar, no conversion needed
    # use the pressLev closest to the bhl to extract the closest geopotential height z_ (phi), closest specific humidity q_ (q) and closest temperature
    geop_blh <- mapply(function(r, c){geop_df[r,c]}, r=1:nrow(geop_df), c=whichCol, SIMPLIFY = T)
    hum_blh <- mapply(function(r, c){hum_df[r,c]}, r=1:nrow(hum_df), c=whichCol, SIMPLIFY = T)
    temp_blh <- mapply(function(r, c){temp_df[r,c]}, r=1:nrow(temp_df), c=whichCol, SIMPLIFY = T)
    # extract missing variables for the calculation. The other variables are surface level once, so we don't have to choose the pressure level of extraction
    mflux <- values(surf_ls[[i]][[grep("ie", names(surf_ls[[i]]))]]) # Instantaneous moisture flux
    sflux <- values(surf_ls[[i]][[grep("ishf", names(surf_ls[[i]]))]]) # Instantaneous surface sensible heat flux
    temp0 <- values(surf_ls[[i]][[grep("t2m", names(surf_ls[[i]]))]]) # Temperature at ground level (2 m), in K
    press0 <- values(surf_ls[[i]][[grep("sp", names(surf_ls[[i]]))]]) # Surface pressure (at ground level), in Pa
    
    # Now we can apply the W* function (convective), specifying the above variables (the constants are within the function)
    # conversions from era5 to the w* equation are taken care of within the function
    w_star <- era5_convective(phi=geop_blh, t_blh=temp_blh, p_blh=pressLev, q_blh=hum_blh, t0=temp0, p0=press0, m_flux=mflux, sh_flux=sflux)
    
    # Put w star values back into the raster of each timestep
    r_wstar <- surf_ls[[i]][[1]]
    values(r_wstar) <- w_star
    names(r_wstar) <- paste0("Wstar_", i)
    
    # Return one layer per timestep
    return(r_wstar)
    
    # Alternatively, we could add W* as extra layer among the surface level variables
    # surf_wstar <- c(surf_ls[[i]], r_wstar)
    # return(surf_wstar)
  })
  
  ## Add wstar as extra variable to the surf level raster stack
  # surf_wstar <- do.call(c, surf_wstar_ls)
  # varnames(surf_wstar) <- c(varnames(surf), "w_star")
  # longnames(surf_wstar) <- c(longnames(surf), "Thermal uplift potential (W*)")
  # # Save the surface variables for each day with the added wstar as a NetCDF file
  # outFile <- gsub("singleLev", "singleLev_Wstar", grep(date, surfFls, value=T))
  # writeCDF(surf_wstar, outFile)
  
  ## Save wstar as separate raster variable, specifying variable names and units
  dailyWstar <- do.call(c, wstar_ls)
  outFile <- paste0(upliftPath,"/download_ERA5_",date,"_Wstar.nc")
  writeCDF(dailyWstar, outFile, varname="w_star", longname="Thermal uplift potential (W*)", unit="m/s", overwrite=TRUE)
  
}, .parallel = T)


#____________________________________________________
# Calculate orographic uplift for each timestep ####
#____________________________________________________

## !!! With a 3 km topography it barely works without crashing. 
# Instead of calculating w_oro for each time step and for the entire world, I simply annotated the final dataset with U-V and with a 1 km DEM, separately, 
# and then computed w_oro directly on the final dataset in the next script!
# With a lower resolution DEM (I tried this successfully with 5 km) or with a computer with more RAM, 
# the following code should calculate w_oro for each hour and each day, 
# matching the spatial scale of the DEM and the temporal scale of the wind data.

#______________________________________________________
# Merge and annotate topography (DEM, slope and aspect)

library(terra)

cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/"
topo_path <- "/home/mscacco/ownCloud/Martina/ENV_DATA/GlobalDEM/"
tmpDir <- tempdir()
lapply(list.files(paste0(topo_path,"globalDEM_voidfill_jonDeFerranti_15sec-450m"), full.names = T), unzip, exdir = tmpDir)
tiles <- list.files(tmpDir, full.names = T, pattern = "tif")

# Merging all at once does not work, let's try to merge portion of the world in equal parts
length(tiles)/8 
#tiles_ls <- split(tiles, cut(seq_along(tiles), breaks=8, labels=FALSE)) # for 30 arc sec, ~1 km res (0.9 at Equator)
tiles_ls <- split(tiles, cut(seq_along(tiles), breaks=2, labels=FALSE)) # for 150 arc sec, ~5 km res (4.6 at Equator)

# create a raster collection and merge of all the tiles
dir.create(paste0(cot_path,"DEM/TOPO_partialMerges_5km"))
lapply(names(tiles_ls), function(i){
  topoColl <- sprc(lapply(tiles_ls[[i]], function(t){
    dem <- rast(t)
    names(dem) <- "dem"
    dem[dem == 0] <- NA # over sea cells get a value of 0. Flatlands over land get a minimum value of 1.
    # dem <- aggregate(dem, fact = 2, fun = "mean") # coarsen resolution to about 1 km (30 arc sec)
    dem <- aggregate(dem, fact = 10, fun = "mean") # coarsen resolution to about 5 km (150 arc sec)
    # slope <- terrain(dem, v="slope", neighbors=8, unit="degrees") # We do this after merging!
    # asp <- terrain(dem, v="aspect", neighbors=8, unit="degrees")
    return(c(dem, slope, asp))
  }))
  merge(topoColl, filename=paste0(cot_path,"DEM/TOPO_partialMerges_5km/worldDEM_5km_part",i,".tif"))
  rm(topoColl); gc()
})
# Merge the partial merges into one large raster:
partMergeColl <- sprc(lapply(list.files(paste0(cot_path,"DEM/TOPO_partialMerges_5km"), full.names = T), rast))
demMerge <- merge(partMergeColl, filename = paste0(cot_path,"DEM/worldDEM_aggMerge_deFerranti_5km.tif"), overwrite=TRUE)
# Slope and aspect need to be calculated after merging, 
# because given the use of neighbouring cells otherwise we get stripes of NA values at the margins of the merged tiles.
demMerge <- rast(paste0(cot_path,"DEM/worldDEM_aggMerge_deFerranti_5km.tif"))
demMerge <- demMerge[["dem"]]
slope <- terrain(demMerge, v="slope", neighbors=8, unit="degrees")
asp <- terrain(demMerge, v="aspect", neighbors=8, unit="degrees")
topoMerge <- c(demMerge, slope, asp)
writeRaster(topoMerge, filename = paste0(cot_path,"DEM/world_DEM_slope_aspect_deFerranti_5km.tif"), overwrite=TRUE)

plot(topoMerge, colNA = "red")


#________________________________________________________________________________
# Add oro uplift calculation as additional variable in the single level variables

setwd("/media/mscacco/Cyprus")
cot_path <- "/home/mscacco/ownCloud/Martina/ProgettiVari/COT/" #desktop
source(paste0(cot_path,"Scripts/COT/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R")) # functions also include necessary libraries
#library(rasterDT)

# This code takes the DEM raster and the surface level wind components (U and V from era5).
# It computes slope and aspect and calculates orographic uplift
# Extract U and V per date
uniqueDates <- sapply(strsplit(list.files("ERA5downloads_singleLev", pattern="nc"), "_"), "[", 3)

# Create a new folder to store the orographic uplift. 
dir.create("ERA5downloads_upliftProxies")
upliftPath <- "ERA5downloads_upliftProxies"

# World topography (dem, slope, aspect):
topoMerge <- rast(paste0(cot_path,"DEM/world_DEM_slope_aspect_deFerranti_5km.tif")) # with the topo at 1 km computation was not possible on this computer
topoDf <- as.data.table(topoMerge, xy=T, na.rm=F) # coords and topo values. Important to not remove the NAs, to make sure all cells get annotated and the values can be put back into the raster template!
nrow(topoDf)==ncell(topoMerge)

# Surface variables file list:
surfFls <- list.files("ERA5downloads_singleLev", full.names = T)

# For each day, we create a raster with hourly W_Oro per pixel (as many layers as the number of hours)
# In this case, given the topography has higher resolution than the era5 surf level variables, the toporgaphy raster will be our template
# Meaning the resulting W_oro hourly rasters will also have a 1 km resolution (instead of ~30 km)
lapply(uniqueDates, function(date){ # per date
  
  # Load rasters for each date
  surf <- rast(grep(date, surfFls, value=T))
  # Split by hour/time-step (both press and surf variables were downloaded for the same hours/time-steps)
  surf_ls <- split(surf, time(surf))
  
  woro_ls <- lapply(seq_along(surf_ls), function(i){ # per hour/time step
    
    # We use the topoRaster as reference, so we use the topo coordinates to extract the U and V values at the locations of the topographic pixels:
    wind_vars <- extract(x = surf_ls[[i]][[grep("u10|v10", names(surf_ls[[i]]))]],
                         y = topoDf[,c("x","y")], method="bilinear")
    
    # # Alternative at lower resolution, using the era5 raster as template, at 30 km:
    # u_wind <- values(surf_ls[[i]][[grep("u10", names(surf_ls[[i]]))]]) # "10 metre U wind component" 
    # v_wind <- values(surf_ls[[i]][[grep("v10", names(surf_ls[[i]]))]]) # "10 metre V wind component"
    # 
    # era5Coords <- crds(surf_ls[[i]])
    # slope <- extract(topoMerge[["slope"]], era5Coords, method="simple")
    # aspect <- extract(topoMerge[["aspect"]], era5Coords, method="simple")

    # Now we can apply the oro uplift function (oroUplift), specifying the above variables
    w_oro <- oroUplift(u=wind_vars[,grep("u10", names(wind_vars))], 
                       v=wind_vars[,grep("v10", names(wind_vars))], 
                       slope=(topoDf$slope * pi / 180), aspect=(topoDf$aspect* pi / 180))
    
    # Put w oro values back into the raster of each timestep
    r_woro <- topoMerge[[1]]
    values(r_woro) <- w_oro
    names(r_woro) <- paste0("woro_", i)
 
    return(r_woro)
    
  })
  
  # Save the orographic uplift for each day as a NetCDF file, specifying variable names and units
  dailyWoro <- do.call(c, woro_ls)
  outFile <- paste0(upliftPath,"/download_ERA5_",date,"_Woro.nc")
  writeCDF(dailyWoro, outFile, varname="w_oro", longname="Orographic uplift potential (Wo)", unit="m/s", overwrite=TRUE)
  
})


