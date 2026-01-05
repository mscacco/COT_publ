library(move2)
library(sf)
library(lubridate)
library(ecmwfr)
library(terra)
library(dplyr)
library(plyr)
library(doMC)
library(data.table)

#________________________________________________
## Functions to request and download era 5 data
#________________________________________________

# Function to be used in following function, to define the area for CDS API
# this function takes either bbox object or a vector of 4 numeric values (in order long min max and lat min max)
area2request <- function(x, ext = 0.01) {
  #check if the input is a bbox
  if(!is(x, "bbox")){
    bounds <- x
  }
  #check whether the input is a vector of length 4 
  if(is(x, "numeric") && length(x) == 4){
    bounds <- data.frame(xmin = x[1], xmax = x[2],  ymin = x[3], ymax = x[4])
  }
  #area: [north, west, south, east]
  x_range <- bounds$xmax - bounds$xmin
  y_range <- bounds$ymax - bounds$ymin
  west <- floor(bounds$xmin - ext * x_range)
  if(west < -180){west <- -180}
  east <- ceiling(bounds$xmax + ext * x_range)
  if(east > 180){east <- 180}
  south <- floor(bounds$ymin - ext * y_range)
  if(south < -90){south <- -90}
  north <- ceiling(bounds$ymax + ext * y_range)
  if(north > 90){north <- 90}
  return(paste(north, west, south, east, sep="/"))
}

# Function to be used in following function, to generate request table from date range
# x is the timestamp column, and unit is for the request (daily request or monthly request)
# For each timestamp, request data for the hour before and after (e.g. for 12:40:00, request data at 12:00:00 and 13:00:00 of that day)
date2request <- function(x, unit = c("day","month")) {
  DateTime <- sort(unique(c(floor_date(x, "hour"), floor_date(x+3600, "hour"))))
  reqTable <- data.frame(year=year(DateTime), month=sprintf("%02d", month(DateTime)), days=sprintf("%02d", day(DateTime)), hour=sprintf("%02d",hour(DateTime)))
  if(unit=="day"){
    #split the table into a list with all the observations with the same year, month, day
    reqTable <- split(reqTable, apply(reqTable[,1:3], 1, paste0, collapse=""))
    #collapse to one line per day
    reqTable <- do.call(rbind, lapply(reqTable, function(x) cbind(x[1,1:3], hours= paste0(x[,4], collapse=" "))))
  }else if(unit == "month"){
    #split the table into a list with all the observations with the same year and month
    reqTable <- split(reqTable, apply(reqTable[,1:2], 1, paste0, collapse=""))
    #collapse to one line per month
    reqTable <- do.call(rbind, lapply(reqTable, function(x) cbind(x[1,1:2], days = paste0(unique(x[,3]), collapse=" "), hours = paste0(unique(x[,4]), collapse=" "))))
  }
  return(reqTable)
}

# Function that uses the above function to compile a request table
# track can be either a move2 object or a data.frame (for a data.frame, the column names of time and coordinates need to be specified unless they fit the default of a Movebank dataframe)
# It is flexible in terms of arranging the request by day or month, and allows the user to select either a predefined area for the entire request ("bbox", useful to then stack the downloaded nc files) or to request different areas per time unit (smaller area downloaded)
# Note that depending on the number of variable requested and pressure levels the monthly unit could create too large files
timeSpace2request <- function(track, 
                              timeUnit = c("day","month"), 
                              area = c("bbox", "perTimeUnit"),
                              timeCol = "timestamp",
                              coordsCol = c("location.long","location.lat")) {
  # Extract time and split data by timeUnit
  if(is(track, "move2")){time <- mt_time(track)}else if(is(track, "data.frame")){time <- track[,timeCol]}
  timeTable <- date2request(x = time, unit = timeUnit)
  # if user selects bbox, the bbox of the entire track is used in the request
  if(area == "bbox"){
    if(is(track, "move2")){x <- st_bbox(track)}else if(is(track, "data.frame")){
      x <- c(min(track[,coordsCol[1]], na.rm=T), 
             max(track[,coordsCol[1]], na.rm=T), 
             min(track[,coordsCol[2]], na.rm=T), 
             max(track[,coordsCol[2]], na.rm=T))
    }
    area <- area2request(x)
    reqTable <- cbind(timeTable, area)
  }else # if user selects "perTimeUnit", the function goes through the timeTable and defines a different bbox per time unit
    if(area == "perTimeUnit"){
    reqTable <- as.data.frame(rbindlist(lapply(1:nrow(timeTable), function(i){
      dayRange <- range(as.numeric(unlist(strsplit(timeTable[i,"days"]," "))))
      minDate <- as.Date(paste(c(timeTable[i,1],timeTable[i,2],dayRange[1]),collapse="-"), tz="UTC")
      maxDate <- as.Date(paste(c(timeTable[i,1],timeTable[i,2],dayRange[2]),collapse="-"), tz="UTC")
      sub <- track[date(time) >= minDate & date(time) <= maxDate,]
      if(is(sub, "move2")){x <- st_bbox(sub)}else if(is(sub, "data.frame")){
        x <- c(min(sub[,coordsCol[1]], na.rm=T), 
               max(sub[,coordsCol[1]], na.rm=T), 
               min(sub[,coordsCol[2]], na.rm=T), 
               max(sub[,coordsCol[2]], na.rm=T))
      }
      area <- area2request(x)
      return(cbind(timeTable[i,],area))
    })))
    }
  # the area to download is added to the time table in both cases, so that dates and area come in the same request table
  return(reqTable)
}


# Function to create a list of requests
# I added a few arguments for additional flexibility in the request
create_request_list <- function(reqTable,
                                datasetName = c("reanalysis-era5-pressure-levels","reanalysis-era5-single-levels"),
                                vars = c("geopotential", "u_component_of_wind", "v_component_of_wind"),
                                levels = c("500", "550", "600", "650", "700", "750", "775", "800", "825",
                                           "850", "875", "900", "925", "950", "975", "1000")) {
  if(datasetName == "reanalysis-era5-single-levels" & length(levels)>0){
    print("Request to download surface (single level) variables, levels are ignored")
    levels <- NULL
    ds <- "singleLev"
  }else{
    levels <- levels
    ds <- "pressLev"}
  reqList <- lapply(seq_len(nrow(reqTable)), function(i) {
    filename <- ifelse(nchar(reqTable$days[i])>2, paste0("download_ERA5_", reqTable$year[i], reqTable$month[i],"_",ds, ".nc"), paste0("download_ERA5_", reqTable$year[i], reqTable$month[i], reqTable$days[i],"_",ds, ".nc"))
    list(
      product_type = "reanalysis",
      variable = vars,
      pressure_level = levels,
      year = reqTable$year[i],
      month = reqTable$month[i],
      day = unlist(strsplit(reqTable$days[i], " ")),
      time = unlist(lapply(strsplit(reqTable$hours[i], " "), paste, ":00", sep = "")),
      area = reqTable$area[i],
      format = "netcdf",
      dataset_short_name = datasetName,
      target = filename
    )
  })
  return(reqList)
}

# Might consider adding a timeout option to avoid problems with curl
send_requests <- function(reqList, localPath, user, batch_size, requestFileName="all_requests.rds") {
  all_requests <- list()
  live_requests <- list()
  dld_requests <- list()
  pause <- 0
  #only submit until all requests are submitted
  while (length(all_requests) < length(reqList)) {
    
    #check whether there is a list of all the requests so far
    #update the status of all requests
    if(file.exists(file.path(localPath, paste0(requestFileName,".rds")))){
      all_requests <- readRDS(file.path(localPath, paste0(requestFileName,".rds")))
      #check for those requests that are listed as completed and should still be available to download and store in live_requests
      live_requests <- try(lapply(all_requests[which(unlist(lapply(all_requests, function(x) x$get_status()))!="deleted")], function(x) x$update_status()), TRUE)
      if(inherits(live_requests, "try-error")){
        live_requests <- list()
        dld_requests <- list()
        unlink(file.path(localPath, paste0(requestFileName,".rds")))
      }else{
        #remove the requests that have already been downloaded and store them in dld_requests
        dld_requests <- all_requests[which(unlist(lapply(all_requests, function(x) x$get_status()))=="deleted")]
      }
    }
    
    # Submit new requests if there are fewer than batch_size active requests
    if(length(live_requests) < batch_size && length(all_requests) < length(reqList)) {
      #submit missing requests
      for(i in 1:(batch_size-length(live_requests))){
        req <- reqList[[length(all_requests) + 1]]
        ncfile <- wf_request(user = user, request = req, transfer = FALSE, 
                             path = localPath, verbose = FALSE) #time_out=
        #store the request in the active_requests list
        live_requests <- append(live_requests, list(ncfile))
        all_requests <- append(all_requests, list(ncfile))
        #adjust the number of submitted requests and active requests
      }
    }
    
    #update the status of the live requests
    live_requests <- lapply(live_requests, function(x) x$update_status())
    #list the ones which are completed
    dldList <- which(unlist(lapply(live_requests, function(x) x$get_status()))=="completed")
    if(length(dldList) == 0){
      message("No requests completed. Pausing for 60 seconds.")
      Sys.sleep(60)
      pause <- pause + 1
      message(paste("Resuming after ", pause, " pauses", ".", sep=""))
    } else {
      #download the completed requests, and add the request calls to the list of downloaded request
      dld_requests <- append(dld_requests, lapply(live_requests[dldList], function(x) {
        x$transfer()
        x$delete()
      }
      ))
      #remove the completed requests from the active requests
      live_requests <- live_requests[-dldList]
    }
    message(paste("Number of active requests:", length(live_requests)))
    #save the list of all requests for resuming broken downloads
    saveRDS(all_requests, file = file.path(localPath, paste0(requestFileName,".rds")))
  }
  
  message(paste("Total requests submitted:", length(all_requests)))
  message(paste("Total files downloaded:", length(dld_requests)))
}


#________________________________________________
## Functions to annotate era 5 data to a dataset
## Interpolating in space, time and press layers
#________________________________________________


# Function to extract relevant data from the bird track and ERA5 file
## Change the function to accept move2 or data.frame as input ##
extract_track_and_era5_data <- function(track, era5_file, 
                                        eventCol = "event.id",
                                        timeCol = "timestamp", 
                                        coordsCol = c("location.long","location.lat"), 
                                        heightCol = "height.above.ellipsoid") {
  # Extract the date from the ERA5 file name
  era5_date <- ymd(as.numeric(strsplit(basename(era5_file), "_")[[1]][3] %>% substr(1, 8)))
  ds <- strsplit(basename(era5_file), "_|\\.")[[1]][4]
  
  # Extract event id (row id), coordinates, times, and heights for the given date
  if(is(track, "move2")){
    event_id <- track$event_id[ymd(paste(year(track$timestamp), sprintf("%02d", month(track$timestamp)), sprintf("%02d", day(track$timestamp)))) == era5_date]
    coords <- as.matrix(st_coordinates(track[ymd(paste(year(track$timestamp), sprintf("%02d", month(track$timestamp)), sprintf("%02d", day(track$timestamp)))) == era5_date, ]))
    times <- mt_time(track[ymd(paste(year(track$timestamp), sprintf("%02d", month(track$timestamp)), sprintf("%02d", day(track$timestamp)))) == era5_date, ])
    heights <- track[ymd(paste(year(track$timestamp), sprintf("%02d", month(track$timestamp)), sprintf("%02d", day(track$timestamp)))) == era5_date, heightCol]
  }
  if(is(track, "data.frame")){
    event_id <- track[ymd(paste(year(track[,timeCol]), sprintf("%02d", month(track[,timeCol])), sprintf("%02d", day(track[,timeCol])))) == era5_date, eventCol]
    coords <- track[ymd(paste(year(track[,timeCol]), sprintf("%02d", month(track[,timeCol])), sprintf("%02d", day(track[,timeCol])))) == era5_date, coordsCol]
    times <- track[ymd(paste(year(track[,timeCol]), sprintf("%02d", month(track[,timeCol])), sprintf("%02d", day(track[,timeCol])))) == era5_date,timeCol]
    heights <- track[ymd(paste(year(track[,timeCol]), sprintf("%02d", month(track[,timeCol])), sprintf("%02d", day(track[,timeCol])))) == era5_date, heightCol]
  }
  # Create a list with the extracted data
  extracted_data <- list(
    era5_file_name = era5_file,
    dataset = ds,
    era5_date = era5_date,
    track_data = data.frame(
      event_id = event_id,
      Long = coords[,1],
      Lat = coords[,2],
      timestamp = times,
      heights = heights
    )
  )
  # Remove rows with missing values
  extracted_data$track_data <- extracted_data$track_data[complete.cases(extracted_data$track_data), ]
  
  return(extracted_data)
}



# Function to annotate the bird track with wind speed and direction data
# Trying to add flexibility in n. of variables to annotate
annotate_era5_data <- function(extracted_data, pathToFolder) { #optional path to save files separately
  # Load the ERA5 raster file
  era5_raster <- rast(extracted_data$era5_file_name)
  nTimeSteps <- length(unique(time(era5_raster)))
  if(nTimeSteps==1){ # give a number to raster variables that only have 1 single time step
    names(era5_raster) <- paste0(names(era5_raster),"_1")
  }
  # if there are locations after 23:00, add to the era5 raster the first hour of the file from the following day
  # if the file doesn't exist, locations after 23:00 will not be annotated (era5 vars == NA)
  if(max(hour(extracted_data$track_data$timestamp)) == 23){ 
    day <- unique(date(extracted_data$track_data$timestamp))
    file_nextDay <- gsub(gsub("-","",day), gsub("-","",day+1), extracted_data$era5_file_name)
    if(!file_nextDay %in% list.files(gsub("/.*","",file_nextDay), full.name=T)){
      warning("Data have locations after 11 pm. Raster file for following date is not available, therefore locations after 11 pm are not annotated.")
    }else{
      era5_nextDay <- rast(file_nextDay)
      era5_nextDay <- era5_nextDay[[hour(time(era5_nextDay)) < 1]]
      # the names of the variables of the first hour of the next day contain a _1 (if there are more than 1 hour) or nothing.
      # To differentiate them from the time steps of the previous day we add a number == to the time steps of previous day + 1.
      idNextDay <- paste0("_",nTimeSteps+1)
      if(any(grepl("_1",names(era5_nextDay)))){
        names(era5_nextDay) <- sub("_1",idNextDay,names(era5_nextDay))
        }else{names(era5_nextDay) <- paste0(names(era5_nextDay),idNextDay)}
      era5_raster <- c(era5_raster, era5_nextDay)
    }
  }
  # Variables to annotate
  vars <- unique(paste0(varnames(era5_raster), "_")) #vars are only one letter; add _ to avoid grepping other colnames including that letter
  vars_longname <- unique(gsub(" ","_",paste0(vars, longnames(era5_raster))))
  
  # Extract coordinates from the input data
  coords <- extracted_data$track_data[, c("Long", "Lat")]
  
  ### Interpolate in space: extract bilinear interpolated data of all variables at the track locations
  interpolated_inSpace <- terra::extract(era5_raster, coords, method = "bilinear")
  interpolated_inSpace <- interpolated_inSpace[, grep(paste(vars, collapse="|"), names(interpolated_inSpace))]
  
  ### Data format and processing changes between the single level and the pressure levels datasets
  ### Surface (single lev) variables have column names of extracted raster with only 2 information: variable_timeStep
  ### Pressure level variables have column names of extracted raster with 3 information: variable_pressLevel_timeStep
  if(grepl("single", extracted_data$dataset)){
    layers_to_process <- unique(vars)
  }else{
    pressLevs <- unique(sapply(strsplit(names(interpolated_inSpace), "_"), "[", 2)) # unique pressure levels
    layers_to_process <- apply(expand.grid(vars, pressLevs), 1, paste, collapse="")
  }
  ### Interpolate in time: calculate trilinear interpolated data, calculate the closest values of the variables based on the track timestamps
  interpolated_inTime <- do.call(cbind, lapply(layers_to_process, function(x) { # for each layer to process and each variable
    layer_id <- grep(x, names(interpolated_inSpace)) 
    tmp <- interpolated_inSpace[, layer_id]
    # each row is one location/observation of a certain variable, each column is a time step
    # for each row (track location)
    # (and, for press lev data, each potential pressure level):
    ret <- data.frame(unlist(lapply(1:nrow(tmp), function(i){ 
      approx(x=unique(time(era5_raster)), # take unique time steps
             y=tmp[i,], # take corresponding values at those time steps (the entire row), for each variable and pressure levels
             xout = extracted_data$track_data$times[i])$y} # interpolate value at exact time of track location
    )))
    names(ret) <- x
    return(ret)
  }))
  # For press lev data, interpolatedInTime will have n.col = var-press level combination. For single level data as many columns as the number of variables.
  
  ### For surface variables nothing more is needed
  if(grepl("single", extracted_data$dataset)){
    # Combine the input data and annotated env data
    names(interpolated_inTime) <- vars_longname # rename more understandably
    annotated_data <- cbind(extracted_data$track_data, interpolated_inTime)
  }
  
  ### For pressure level variables we interpolate in height: interpolate the value of the variables at the closest pressure level, depending on the track height at each location
  ### Each level has a constant pressure, but the corresponding height in m vary depending on the Geopotential (gravitational potential energy of a unit mass at a particular location m2/s2, relative to mean sea level). 
  ### We calculate it below by dividing the geopotential by the Earth's gravitational acceleration, g (=9.80665 m s-2).
  if(grepl("press", extracted_data$dataset)){
    interpolated_inHeight <- do.call(cbind, lapply(vars[-grep("z_",vars)], function(v){
      interpolatedVar <- data.frame(unlist(lapply(1:nrow(interpolated_inTime), function(x){
        approx(interpolated_inTime[x, grep("z_", names(interpolated_inTime))] / 9.80665, #geopotential height
               interpolated_inTime[x, grep(v, names(interpolated_inTime))], #var
               xout = extracted_data$track_data$heights[x])$y # one interpolated value per location
      })))# the output is one interpolated value for each track location/time/height
      names(interpolatedVar) <- vars_longname[grep(v, vars)]
      return(interpolatedVar)
    })) # for each var
    
    # Combine the input data and annotated env data
    annotated_data <- cbind(extracted_data$track_data, interpolated_inHeight)
  }
  
  if(!missing(pathToFolder)){
    saveRDS(annotated_data, file = file.path(pathToFolder, paste0("annotatedData_",extracted_data$era5_date, "_", extracted_data$dataset, ".rds")))
  }
  return(annotated_data)
}


#________________________________________________________________
## Functions to calculate thermal and orographic uplift velocity
#________________________________________________________________

#__________________________
# Thermal uplift (w star):

# a function to calculate w* from Hester Bronnvik 2024 (https://www.cell.com/current-biology/fulltext/S0960-9822(24)00392-0)
#I ended up including this bit by bit in the final dplyr pipeline, as running the function as is uses too much RAM
era5_convective <- function(phi, t_blh, p_blh, q_blh, t0, p0, m_flux, sh_flux) {
  
  # Where:
  # phi is the geopotential at (closest to) the boundary layer height
  # q_blh is the specific humidity at the boundary layer height (kg/kg)
  # p_blh is the pressure at the boundary layer height, that is the closest pressure level (in hPa or mbar)
  # t_blh is the temperature below the boundary layer height (in K)
  # t0 is the temperature at ground level (in K, to be converted to Celsius)
  # p0 is the pressure at ground level (in Pa, to be converted to mbar)
  # sh_flux is the surface sensible heat flux (negative upwards, to be reversed)
  # m_flux is the moisture flux (negative upwards, to be reversed)
  
  # Conversions from ERA5 standards units to the units for the calculation:
  t0 <- t0 - 273.15 # temperature at ground level, from K to Celsius
  p0 <- p0 / 100 # pressure at ground level, from Pa to mbar
  sh_flux <- sh_flux * -1 # reverse the sign. ECMWF upward fluxes are negative
  m_flux <- m_flux * -1
  
  # Define constants:
  k <- 0.2854 # Poisson constant
  c_p <- 1012 # the isobaric mass heat capacity of air at common conditions
  p <- 1.225 # the density of air at sea level and 273 K
  
  # Calculate W*:
  Theta_k_z <- t_blh * ((p0 / p_blh) ^ k)
  Thetav_k_z <- Theta_k_z * (1+(0.61 * q_blh))
  wT <- sh_flux / c_p / p 
  wq <- m_flux / p 
  
  wthetav <- wT + 0.61 * t0 * wq
  
  w_star <- .CubeRoot(phi * (wthetav / Thetav_k_z))
  return(as.numeric(w_star))
  
}

# a function to take the cube root
.CubeRoot <- function(x){
  sign(x)*abs(x)^(1/3)
}

#____________________________
# Orographic uplift (w oro):

# ## Orographic uplift test
# asp <- seq(0,360,45) * pi / 180
# b <- seq(0,360,45) * pi / 180
# slop <- seq(0,45,5) * pi / 180
# ws <- seq(0,40,10)
# test <- expand.grid(asp, b, slop, ws)
# names(test) <- c("aspect","beta","slope","windSpeed")
# summary(sin(test$slope))
# test$Ca <- sin(test$slope) * cos(test$aspect-test$beta)
# test$wOro <- test$windSpeed * test$Ca
# summary(test[,c("Ca","wOro")])
# test[test$Ca < 0,]


## Aspect angles conversions:

angles_to_cartesian <- function(x) { # x is a raster of angles
  sin_values <- sin(x * pi / 180)
  cos_values <- cos(x * pi / 180)
  cartesianRast <- c(sin_values, cos_values)
  return(cartesianRast)
}

cartesian_to_angles <- function(x_sin, x_cos) { # sin and cosin from which we want to obtain angles
  angles <- atan2(x_sin, x_cos) * 180 / pi
  return((angles + 360) %% 360) # Ensure angles are within 0-360 degrees
}

## w_oro calculation:

oroUplift <- function(u, v, slope, aspect) {
  
  # Where:
  # u and v are the two vector components of wind at 10 m above the surface (era5 single level variables)
  # slope is the angle of the slope in radians (theta)
  # aspect is the aspect of the slope in radians between 0 and 2pi (alpha)
  
  # Calculate horizontal wind speed (ws) and wind direction FROM (beta):
  ws <- sqrt(u^2 + v^2) # wind speed at 10 m height, in m/s
  beta <- atan2(u, v) + pi # wind dir at 10 m, blowing FROM, in rad between 0 and 2pi

  # Calculate Ca, the lifting coefficient (if the slope value is 0, the whole Ca coefficient will be 0, so flat slopes have no oro uplift): Ca <- sin(theta) * cos(alpha - beta)
  Ca <- sin(slope) * cos(aspect - beta) 
  # Calculate the orographic uplift according to Brandes et al (2004) and Bohrer et al (2012)
  oroUpl <- ws * Ca
  
  return(oroUpl)
}

windSpeed <- function(u, v) {
  return(sqrt(u^2 + v^2))
}
# With this function we calculate the wind direction (blowing FROM) in degrees from 0 to 360 from the NORTH (y axis)
# For the last step of this calculation see discussion at: https://stackoverflow.com/questions/21484558/how-to-calculate-wind-direction-from-u-and-v-wind-components-in-r
## N.B. if not transformed, ATAN2(U,V) RETURNS THE WIND DIRECTION IN RADIANS FROM -PI to +PI
# For a slope facing eastward (alpha = 90deg), the deflection is maximal when the wind comes from the East and minimal (negative) when the wind comes from the West. 
# This means that the wind direction has to be the meteorological wind (the direction FROM where the winds blows). 
# This way, alpha - beta = 0 and the cosine is maximal (=1).
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
windDir_from_deg <- function(u, v) { 
  WD.deg <- atan2(u, v) * 180/pi  #from rad to deg as +/-180
  return((WD.deg + 180) %% 360) #reverse from blowing "to" to blowing FROM, as 0-360. The %% 360 makes sure the results stays within 0-360
}

windDir_from_rad <- function(u, v) { 
  return(atan2(u, v) + pi) # adding pi rotates the angle in radians and at the same time shift them into the 0-2pi range
}
#calculate wind support (if negative head wind, if positive tail wind)
# input: u and v are the wind components you download, dg is the track direction and has to range from 0 to 360
windSupport <- function(u, v, dg){ 
  if(any(dg[!is.na(dg)] < 0)){stop("The track direction (Dg) has negative values, it has to range from 0 to 360!")}
  dg_rad <- dg/180*pi                        #transform Dg in radians from 0 to 2*pi
  wd <- atan2(u, v)                          #wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)    #wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                    #beta is the abs of wd-dg (they are both from 0 to 2*pi)
  return(cos(beta) * sqrt(u * u + v * v))    #ws = cos(beta)*Vw
}
# apply the function as:
# ws <- wind.support(data$u, data$v, data$dg)

#calculate cross wind, absolute value so no different if from one side or the other side of the animal
#input: u and v are the wind components you download, dg is the track direction and has to range from 0 to 360

crossWind <- function(u, v, dg){
  if(any(dg[!is.na(dg)] < 0)){stop("The track direction (Dg) has negative values, it has to range from 0 to 360!")}
  dg_rad <- dg/180*pi                           #transform Dg in radians from 0 to 2*pi
  wd <- atan2(u, v)                             #wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)       #wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                       #beta is the abs of wd-dg (they are both from 0 to 2*pi)
  return(abs(sin(beta) * sqrt(u * u + v * v)))  #wc = |sin(beta)*Vw|
}
# apply the function as:
# cw <- cross.wind(data$u, data$v, data$dg)

#calculate airspeed, you need cross wind, wind support and track ground speed (Vg)
#Va = sqrt((Vg-Ws)^2 + Wc^2)

airspeed <- function(Vg, Ws, Cw) {
  return(sqrt((Vg - Ws)^2 + (Cw)^2))
}
