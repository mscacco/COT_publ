#__________________________________________________________
# Bind all studies/datasets together in one data.frame ####
#__________________________________________________________

# This combined data.frame will have one entry per gps fix, only for commuting segments

# Set as wd the path where you stored the content downloaded from the Edmond repository
# which should also be the parent folder where you stored all intermediate results of previous scripts.
setwd("...")


#__________________________
# Standardise the datasets of all studies to bind them together

library(scales)
library(data.table)

fls <- list.files("DataProcessed", "_dfGpsAcc_allSegmentsID_thresholdClass", full.names=T)

# Check the column names of all studies and see if they match to bind the studies
allColLs <- lapply(fls, function(f){
  load(f)
  return(names(df_allSegm))
})
allColNames <- do.call(c, allColLs)
# The column that are contains in only some studies are not important for us so we can remove them
table(allColNames)[table(allColNames)==length(fls)]
table(allColNames)[table(allColNames)<length(fls)]
# And keep only those contained in all studies + other important columns
colsToKeep <- c("timestamp", "acc_sampl_freq_per_axis","meanVedba_Gs","cumVedba_Gs", names(table(allColNames)[table(allColNames)==length(fls)]))
colsToAdd <- c("acc_closest_timestamp","acc_event_id", "diff_acc_time_s","acc_axes", "acc_Naxes", "cumVedba", "flapping",
               "deployment.id","individual.id","tag.id","individual.tag.id","study.id","tag.local.identifier")
heightCols <- c("height","height.above.ellipsoid","height.above.msl","height.msl")
sampFreqCols <- c("acc_sampl_freq","acc_sampl_freq_per_axis")
#fls[!sapply(sapply(allColLs, "%in%", sampFreqCols), any)]

ellp <- fls[sapply(allColLs, function(c){"height.above.ellipsoid" %in% c})]
raw <- fls[sapply(allColLs, function(c){"height.raw" %in% c})] 
raw %in% ellp #fls with height.raw also have height.above.ellipsoid, so we can ignore the height.raw column
fls[sapply(allColLs, function(c){"height.msl" %in% c})] #allison's data
fls[sapply(allColLs, function(c){"height" %in% c})] #fls with height are from DDiaries, Ornitela and Technosmart tags, these are height msl


## Now we can bind all studies, keeping only commuting (flight) segments
# there is a warning of NA coercion, to know where it comes from, turn warnings into errors:
#If warn is 2 all warnings are turned into errors. If warn is negative all warnings are ignored. 
#If warn is 0 (the default) warnings are stored until the top–level function returns.
options(warn=0) #options(warn=2)
allStudies_ls <- lapply(fls, function(f)try({
  print(f)
  load(f) #object df_allSegm
  studySpeciesId <- paste(strsplit(f, "_")[[1]][2:3], collapse="_")
  deviceType <- strsplit(f, "_")[[1]][4]
  if(all(sampFreqCols %in% names(df_allSegm)==F)){
    df_allSegm$acc_sampl_freq_per_axis <- df_allSegm$n_samples_per_axis / df_allSegm$acc_burst_duration_s
  }else if("acc_sampl_freq" %in% names(df_allSegm)){
    names(df_allSegm)[names(df_allSegm)=="acc_sampl_freq"] <- "acc_sampl_freq_per_axis"
  }
  # time column
  if("timestamps" %in% names(df_allSegm)){
    names(df_allSegm)[names(df_allSegm)=="timestamps"] <- "timestamp"
  }
  # Keep only common columns
  if(any(colsToAdd %in% names(df_allSegm)==F)){
    colsSub <- colsToAdd[!colsToAdd %in% names(df_allSegm)]
    df_allSegm[,colsSub] <- NA
  }
  # Add vedba columns in Gs where missing
  if(!"meanVedba_Gs" %in% names(df_allSegm)){ #all non-eobs tags are measured in Gs, so we create a new column and duplicate the existing one
    df_allSegm$meanVedba_Gs <- df_allSegm$meanVedba
    df_allSegm$cumVedba_Gs <- df_allSegm$cumVedba
  }
  # Name height as either above ellipsoid or above sea and create a generic height column that contains either information
  df_allSegm$height.raw <- NULL
  heightCol <- grep("height", names(df_allSegm), value=T)
  if(length(heightCol)==0){
    df_allSegm$height_gener <- NA
  }
  if(length(heightCol)>0){
    if(length(heightCol)==2){heightCol <- heightCol[which.min(c(length(which(is.na(df_allSegm[,heightCol[1]]))),length(which(is.na(df_allSegm[,heightCol[2]])))))]}
    df_allSegm$height_gener <- df_allSegm[,heightCol]
    if(heightCol %in% c("height","height.msl")){df_allSegm$height.above.msl <- df_allSegm[,heightCol]}
  }
  if(!"height.above.ellipsoid" %in% names(df_allSegm)){df_allSegm$height.above.ellipsoid <- NA}
  if(!"height.above.msl" %in% names(df_allSegm)){df_allSegm$height.above.msl <- NA}
  df_allSegm <- df_allSegm[,c(colsToKeep, colsToAdd, "height.above.ellipsoid", "height.above.msl", "height_gener")]
  names(df_allSegm)[names(df_allSegm) %in% c("coords.x1","coords.x2")] <- c("location.long","location.lat")
  # Coerce columns of all studies to the same classes (all characters)
  df_allSegm$diff_acc_time_s <- as.numeric(df_allSegm$diff_acc_time_s)
  df_allSegm[c("acc_closest_timestamp","acc_event_id","individual.tag.id","tag.local.identifier","deployment.id","individual.id","tag.id","study.id")] <- 
    sapply(df_allSegm[c("acc_closest_timestamp","acc_event_id","individual.tag.id","tag.local.identifier","deployment.id","individual.id","tag.id","study.id")],as.character)
  # Keep only flight segments
  df_flightSegm <- df_allSegm[!is.na(df_allSegm$track_flight_id),]
  if(nrow(df_flightSegm)>0){
    # Coerce factor columns to character before binding the dataframes (otherwise factor columns will have non-matching levels)
    factorCols <- names(df_flightSegm)[which(sapply(df_flightSegm, class)=="factor")]
    df_flightSegm[,factorCols] <- sapply(df_flightSegm[,factorCols], as.character)
    # Save histogram of the mean Vedba
    png(paste0("Plots/finalFlightSegments_thresholdClassification/selectedCommutingSegments_meanVedba/study_",studySpeciesId,"_",deviceType,"_finalSelectedSegments_MeanVedbaHistogram.png"),
        width=8,height=6,units="in",res=300)
    hist(df_flightSegm$meanVedba_Gs, breaks = "FD", col="grey", xlab="mean VeDBA (Gs)", main="")
    dev.off()
    # return the subset dataframe of flight segments to the list
    return(df_flightSegm)
  }else if(nrow(df_flightSegm)==0){warning("No observations were classified as commuting flight, study gets excluded.")}
}))

# check errors
is.error <- function(x) inherits(x, "try-error")
errors <- vapply(allStudies_ls, is.error, logical(1))
allStudies_ls[errors]
fls[errors]
table(sapply(allStudies_ls, class))
allStudies_ls[sapply(allStudies_ls, class)!="data.frame"]
allStudies_ls <- allStudies_ls[sapply(allStudies_ls, class)=="data.frame"]
allStudies <- as.data.frame(rbindlist(allStudies_ls, use.names=T))

# make sure coordinates and timestamp and sampling frequency are not missing
c("timestamp","location.long","location.lat","meanVedba_Gs","acc_sampl_freq_per_axis","acc_Naxes") %in% names(allStudies)

# Some other checks
length(unique(allStudies$individual.taxon.canonical.name))
length(unique(allStudies$study.id))
length(unique(allStudies$study.name))
table(allStudies$deviceType)

# add a column that separates movebank from external studies
allStudies$dataSource <- "Movebank"
allStudies$dataSource[is.na(allStudies$study.id)] <- "External"
unique(allStudies$study.name[allStudies$dataSource == "External"])

# Call all the Technosmart tags with the same name
table(allStudies$deviceType)
allStudies$deviceType[allStudies$deviceType == "Technosmart_axytreck"] <- "Technosmart"
table(allStudies$deviceType)

# Save final binded dataframe
save(allStudies, file="DataFinalSummary/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024.RData")

#______________________________________________________________________
# Check and remove duplicated individuals across Movebank studies ####

# In Movebank it is possible that a same individual has been stored twice in different studies, and potentially with different names
# So it is important to check both individual id and tag id
# if an individual/tag combination is duplicated, we make also check for the species (to make sure it is in fact the same individual) and the time interval
# same individual in different studies with non-overlapping times are kept

load("DataFinalSummary/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024.RData") #object allStudies

# Check for duplicated individuals between MOVEBANK studies
# remove studies that are not in movebank (all those that do not have a study.id)
allStudies$individual.study.id <- paste(allStudies$individual.local.identifier,allStudies$study.id, sep="_")
unique(allStudies$study.name[is.na(allStudies$study.id)]) # they are all non-movebank studies so not affected by duplicates
allStudies_movebank <- allStudies[which(allStudies$dataSource == "Movebank"),]
# we make a dataframe with one entry per individuals per study, including the time range covered by each individuals and the number of entries
dupIndStud <- duplicated(allStudies_movebank[,c("individual.local.identifier","study.id")])
unique_IndStudies <- allStudies_movebank[!dupIndStud, c("individual.local.identifier","tag.local.identifier","tag.id","study.id","study.name","individual.study.id","individual.taxon.canonical.name")]
# length(unique(allStudies_movebank$study.name))
# length(unique(unique_IndStudies$study.name))
rownames(unique_IndStudies) <- 1:nrow(unique_IndStudies)
table(duplicated(unique_IndStudies[,c("individual.local.identifier")]))
unique_IndStudies[duplicated(unique_IndStudies[,c("individual.local.identifier")]),]
table(duplicated(unique_IndStudies[,c("individual.study.id")]))
# there are no duplicates in individual.study.id, so we can group by it
indTimeRange <- data.frame(aggregate(timestamp~individual.study.id, data=allStudies_movebank, FUN=min),
                           timestampEnd=aggregate(timestamp~individual.study.id, data=allStudies_movebank, FUN=max)[,2])
Nsegments <- as.data.frame(table(allStudies_movebank$individual.study.id))
unique_IndStudies <- merge(unique_IndStudies, indTimeRange, by="individual.study.id", all.x = T)
unique_IndStudies <- merge(unique_IndStudies, Nsegments, by.x="individual.study.id", by.y="Var1", all.x = T)
names(unique_IndStudies)[8:10] <- c("timestampStart","timestampEnd","Nsegments")
unique_IndStudies$timestampStart <- as.Date(trunc(unique_IndStudies$timestampStart, "days"))
unique_IndStudies$timestampEnd <- as.Date(trunc(unique_IndStudies$timestampEnd, "days"))
# now we can check duplicated individuals across studies, there are no duplicated individual names
# unique_IndStudies[unique_IndStudies$individual.taxon.canonical.name=="Calonectris diomedea",]
###---There is no duplicates based on the individual name, but there are duplicates based on the tag id:---###
unique_IndStudies <- unique_IndStudies[!is.na(unique_IndStudies$tag.local.identifier),]
table(duplicated(unique_IndStudies[,c("tag.local.identifier","individual.taxon.canonical.name")]))
unique_IndStudies$study.name[duplicated(unique_IndStudies[,c("tag.local.identifier","individual.taxon.canonical.name")])]
# create a new variable that associates the tag.id with the species
unique_IndStudies$tag.species.identifier <- paste(unique_IndStudies$individual.taxon.canonical.name, unique_IndStudies$tag.local.identifier, sep="_")
dupTags <- unique_IndStudies[duplicated(unique_IndStudies[,c("tag.species.identifier")]),"tag.species.identifier"]
dupTags <- unique(dupTags)
# Extract all rows from the unique identifier data that have duplicated tag.local.identifier
dataDuplTags <- unique_IndStudies[which(unique_IndStudies$tag.species.identifier %in% dupTags),]
# Split by tag ID so we will have in each element duplicated of the same information
duplTags_ls <- split(dataDuplTags, dataDuplTags$tag.species.identifier)
# Each element of the list should only belong to one species
table(sapply(duplTags_ls, function(tag) length(unique(tag$individual.taxon.canonical.name))))
# Keep only duplicated that occur across studies (if within the same study it is probably the same tag on different individuals)
duplTags_ls <- lapply(duplTags_ls, function(tag){
  if(length(unique(tag$study.id))>1){
    return(tag)
  }
})
duplTags_ls <- duplTags_ls[!sapply(duplTags_ls, is.null)]
# Now in each of the elements left there are some duplicated tags, we need to choose which to keep
# Check how many duplicates there are
table(sapply(duplTags_ls, nrow))
# Some have more than 2 duplicates, so we write a function that checks and chooses between multiple duplicates
library(DescTools) # for function overlap
duplTags_ls <- lapply(duplTags_ls, function(tag){
  tag$keep <- TRUE
  # set to keep=true all rows and modify to false the individuals that we want to exclude 
  # (in case of overlap in time, we exclude the individuals with less segments)
  # in case of no overlap we keep them all
  # we write a general code that accounts for multiple duplicates:
  rowCombs <- expand.grid(row.names(tag), row.names(tag), stringsAsFactors = F)
  rowCombs <- rowCombs[rowCombs$Var1!=rowCombs$Var2,]
  rowCombs$overlap <- FALSE
  # Check if each combination of individuals overlaps in time
  for(i in 1:nrow(rowCombs)){
    if(as.numeric(tag[rowCombs[i,"Var1"],c("timestampStart","timestampEnd")]) %overlaps% 
       as.numeric(tag[rowCombs[i,"Var2"],c("timestampStart","timestampEnd")]) == TRUE){
      rowCombs[i,"overlap"] <- TRUE
    }}
  # For those which overlap exclude (keep=F) the one with less segments
  dupRows <- rowCombs$Var1[rowCombs$overlap==T]
  toExclude <- dupRows[which.min(tag[dupRows,"Nsegments"])]
  tag[row.names(tag)==toExclude, "keep"] <- FALSE
  # with only two would have been easier:
  # if(nrow(tag)==2){
  #   Overlap(as.numeric(tag[1,c("timestampStart","timestampEnd")]), as.numeric(tag[2,c("timestampStart","timestampEnd")]))
  #   if(as.numeric(tag[1,c("timestampStart","timestampEnd")]) %overlaps% 
  #      as.numeric(tag[2,c("timestampStart","timestampEnd")]) == TRUE){
  #     tag$keep[which.max(tag$Nsegments)] <- TRUE #if both have same number of segments it will automatically take the first occurrence
  # }
  return(tag)
})
duplTags_df <- rbindlist(duplTags_ls)

# Individuals to Remove from the original dataset (the one loaded at the beginning of the script, containing both Movebank and non Movebank studies):
# (first check that individual.study.id is unique)
length(unique(duplTags_df$individual.study.id)) == nrow(duplTags_df)
individualsToRemove <- duplTags_df$individual.study.id[duplTags_df$keep == F]

allStudies_noDups <- allStudies[!allStudies$individual.study.id %in% individualsToRemove,]

# for external studies for which individual.tag.id column does not exist, create it and copy the content of individual.local.identifier
allStudies_noDups$individual.tag.id[which(allStudies_noDups$dataSource == "External")] <- allStudies_noDups$individual.local.identifier[which(allStudies_noDups$dataSource == "External")]

# some numeric checks
length(unique(allStudies_noDups$individual.taxon.canonical.name))
length(unique(allStudies_noDups$individual.tag.id))
length(unique(allStudies_noDups$study.name))
length(unique(allStudies_noDups$study.name[allStudies_noDups$dataSource=="External"]))
table(allStudies_noDups$dataSource)

#_______________________________________________________
# Quick sanity checks and data standards summaries ####

# In madebytheo tags, all bursts are 20 Hz, and some bursts only contain 0 for all values in the burst
# This only happens for one species in one study, below you can see the data exploration.
# Looking at the distribution of the data it looks like it's a sensor error, so we remove these data.
nrow(allStudies_noDups[which(allStudies_noDups$meanVedba==0),])
table(allStudies_noDups$acc_sampl_freq_per_axis[which(allStudies_noDups$meanVedba==0)])
table(allStudies_noDups$deviceType[which(allStudies_noDups$meanVedba==0)])
unique(allStudies_noDups[which(allStudies_noDups$meanVedba==0),c("study.name","study.id","individual.taxon.canonical.name")])
test <- allStudies_noDups[allStudies_noDups$study.id==180156318 & allStudies_noDups$individual.taxon.canonical.name == "Anser albifrons",]
hist(test$meanVedba_Gs, breaks="FD")
range(test$meanVedba_Gs, na.rm=T)
length(test$meanVedba_Gs[which(test$meanVedba_Gs == min(test$meanVedba_Gs, na.rm=T))])
length(test$meanVedba_Gs[which(test$meanVedba_Gs == max(test$meanVedba_Gs, na.rm=T))])
minmax <- range(test$meanVedba_Gs, na.rm=T)
toRemove <- which(allStudies_noDups$study.id == 180156318 & 
                    allStudies_noDups$individual.taxon.canonical.name == "Anser albifrons" &
                    allStudies_noDups$meanVedba_Gs %in% minmax)
allStudies_noDups <- allStudies_noDups[-toRemove,]

# Remove missing values in Vedba (only 25 observations)
table(complete.cases(allStudies_noDups$meanVedba_Gs))
allStudies_noDups <- allStudies_noDups[complete.cases(allStudies_noDups$meanVedba_Gs),]

#____________
# Save dataset without duplicates. This dataset is stored in the Edmond repository
saveRDS(allStudies_noDups, file="BiologgingData/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024_noDupl.rds")

#__________
# Summaries

# Some descriptive stats on the acc sampling schedules we have
table(allStudies_noDups$acc_axes)
table(allStudies_noDups$acc_sampl_freq_per_axis)
summary(allStudies_noDups$acc_sampl_freq_per_axis)
summary(allStudies_noDups$acc_burst_duration_s)
summary(allStudies_noDups$n_samples_per_axis)

overView_ACCsamplFreq <- table(allStudies_noDups$individual.taxon.canonical.name, 
                               allStudies_noDups$acc_sampl_freq_per_axis)

length(unique(allStudies_noDups$individual.taxon.canonical.name))
unique(allStudies_noDups$study.name[allStudies_noDups$dataSource=="External"])

# summary of the time difference between gps and acc
summary(allStudies_noDups$diff_acc_time_s)

# Check the distribution of meanVedba per device type
table(allStudies_noDups$deviceType)
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="eobs"], breaks="FD", main="eobs")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="Ornitela"], breaks="FD", main="Ornitela")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="Technosmart"], breaks="FD", main="Technosmart")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="DailyDiary"], breaks="FD", main="Daily Diary")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="milsar"], breaks="FD", main="Milsar")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="madebytheo"], breaks="FD", main="madebytheo")
hist(allStudies_noDups$meanVedba_Gs, breaks="FD", main="All devices")


