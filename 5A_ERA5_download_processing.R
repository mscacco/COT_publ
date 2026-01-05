
#________________________
# Download ERA5 data ####
#________________________

## Download environmental variables (ERA5 data)
## Script and functions modified from: https://github.com/kamransafi/ERA5_move2
## We will download and annotate U_wind, V_wind, and all variables needed for the calculation of oroUplift and thermUplift
## To run this script the dataset resulting from the script "3B.VEDBAclassification_MixR_createFinalDf.R" is needed


setwd("...") # to store env data
dir.create("ERA5downloads")
localPath <- "ERA5downloads"

cot_path <- "..." # local folder with data and scripts
source(paste0(cot_path,"Scripts/COT_publ/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R")) # load own function library
finalDf <- readRDS(paste0(cot_path,"DataFinalSummary/FinalDf_perPoint_VedbaGs_flappingProbs.rds"))

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
    
  })

  ## Save wstar as separate raster variable, specifying variable names and units
  dailyWstar <- do.call(c, wstar_ls)
  outFile <- paste0(upliftPath,"/download_ERA5_",date,"_Wstar.nc")
  writeCDF(dailyWstar, outFile, varname="w_star", longname="Thermal uplift potential (W*)", unit="m/s", overwrite=TRUE)
  
}, .parallel = T)




