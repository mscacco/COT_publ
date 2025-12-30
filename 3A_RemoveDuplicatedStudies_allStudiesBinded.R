#__________________________________________________________
# Bind all studies/datasets together in one data.frame ####
#__________________________________________________________
# This will have one entry per gps fix, only for commuting segments

setwd("...")

#_________
# Add this small step at the beginning, applied to ALL files, just to add a column containing the number of ACC axes (we forgot to add it in step 1A)

fls <- list.files("DataProcessed", "_dfGpsAcc_allSegmentsID_thresholdClass", full.names=T)
fls_acc <- list.files("DataDownloaded", "onlyAcc", full.names=T)

err <- lapply(fls, function(f)try({
  print(f)
  studyId <- strsplit(f, "_")[[1]][2]
  load(f) #object df_allSegm
  if(grepl("Flavio|David|Allison|Emily|Paolo", f)){ #These external datasets (processed in scripts 1AB, dailydiaries and technosmart tags) all collected ACC on 3 axes
    df_allSegm$acc_axes <- "XYZ"
    df_allSegm$acc_Naxes <- 3
  }
  if(!"acc_axes" %in% names(df_allSegm)){
    load(grep(studyId, fls_acc, value=T)) #data.frame, object accDf
    accDf <- as.data.frame(accDf)
    axesCol <- grep("acceleration.axes|acceleration_axes", names(accDf), value=T)
    if(length(axesCol)==2){
      axesCol <- axesCol[which.min(c(length(which(is.na(accDf[,axesCol[1]]))), length(which(is.na(accDf[,axesCol[2]])))))] #compare the number of NAs and exclude the one with more
    } # create two NA axes column (default)
    df_allSegm$acc_axes <- NA
    df_allSegm$acc_Naxes <- NA
    if("acc_event_id" %in% names(df_allSegm) | length(axesCol)>0){ #some datasets don't have an acc_event_id column or an axes column
      if(length(unique(df_allSegm$acc_event_id[!is.na(df_allSegm$acc_event_id)]))!=nrow(df_allSegm[!is.na(df_allSegm$acc_event_id),]) | 
         length(unique(accDf$event_id))!=nrow(accDf)){
        if(length(unique(accDf[,axesCol]))==1){
          df_allSegm$acc_axes <- as.character(unique(accDf[,axesCol]))
          df_allSegm$acc_Naxes <- nchar(as.character(unique(accDf[,axesCol])))
        }
      }else{
        df_allSegm <- merge(df_allSegm, accDf[,c("event_id", axesCol)], by.x="acc_event_id", by.y="event_id", all.x=T)
        names(df_allSegm)[ncol(df_allSegm)] <- "acc_axes"
        df_allSegm$acc_axes <- as.character(df_allSegm$acc_axes)
        df_allSegm$acc_Naxes <- nchar(df_allSegm$acc_axes)
      }
    }
    save(df_allSegm, file=f) #overwrite existing file with added column
  }
}))

# Check errors
is.error <- function(x) inherits(x, "try-error")
errors <- vapply(err, is.error, logical(1))
fls <- fls[51:length(fls)][errors]


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
# Check and remove duplicated individuals across movebank studies ####

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

save(allStudies_noDups, file="DataFinalSummary/allStudies_allTags_allFlightSegments_binded_noDuplicatedIndividuals_birdsBats_thresholdClass_transfGs_March2024.RData")

#_______________________________________________________
# Quick sanity checks and data standards summaries ####
#_______________________________________________________

load("DataFinalSummary/allStudies_allTags_allFlightSegments_binded_noDuplicatedIndividuals_birdsBats_thresholdClass_transfGs_March2024.RData") #object allStudies_noDups

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

# Replace meanVedba plot after removing vedba 0s:
anser <- allStudies_noDups[allStudies_noDups$study.id == 180156318 & allStudies_noDups$individual.taxon.canonical.name == "Anser albifrons",]
png("Plots/finalFlightSegments_thresholdClassification/selectedCommutingSegments_meanVedba/study_180156318_Anser-albifrons_madebytheo_finalSelectedSegments_MeanVedbaHistogram.png",
    width=8,height=6,units="in",res=300)
hist(anser$meanVedba_Gs, breaks = "FD", col="grey", xlab="mean VeDBA (Gs)", main="")
dev.off()

# Remove missing values in Vedba (only 25 observations)
table(complete.cases(allStudies_noDups$meanVedba_Gs))
allStudies_noDups <- allStudies_noDups[complete.cases(allStudies_noDups$meanVedba_Gs),]

# Re-save the dataset without these 430 observations that had meanVedba at their min and max (probably because at the extremes of what the sensor could collect)
save(allStudies_noDups, file="DataFinalSummary/allStudies_allTags_allFlightSegments_binded_birdsBats_thresholdClass_transfGs_March2024_noDupl.RData")

# Now make a table with the final sampling frequencies
table(allStudies_noDups$acc_sampl_freq_per_axis)
summary(allStudies_noDups$acc_sampl_freq_per_axis)
# Only the oilbirds have 1 Hz, but for now we keep them
table(allStudies_noDups$deviceType[which(allStudies_noDups$acc_sampl_freq_per_axis==1)])
table(allStudies_noDups$study.name[which(allStudies_noDups$acc_sampl_freq_per_axis==1)])

overView_ACCsamplFreq <- table(allStudies_noDups$individual.taxon.canonical.name, 
                               allStudies_noDups$acc_sampl_freq_per_axis)
write.csv(overView_ACCsamplFreq, "DataFinalSummary/overview_ACCsamplingFreq_perSpecies_thresholdClass.csv", row.names = F)

length(unique(allStudies_noDups$individual.taxon.canonical.name))
unique(allStudies_noDups$study.name[allStudies_noDups$dataSource=="External"])
length(unique(allStudies_noDups$individual.taxon.canonical.name[allStudies_noDups$acc_sampl_freq_per_axis>8]))
nrow(allStudies_noDups[which(allStudies_noDups$meanVedba==0 & allStudies_noDups$acc_sampl_freq_per_axis>8),])

# Filter out observations with ACC data with very low sampling frequency (less than 8 points in a second)
#allStudies_noDups <- allStudies_noDups[which(allStudies_noDups$acc_sampl_freq_per_axis>8),]

# time difference from gps goes from 0 to 5 min, just to keep in mind
summary(allStudies_noDups$diff_acc_time_s)

#plot(meanVedba_Gs~acc_sampl_freq_per_axis, data=allStudies_noDups)

# Some descriptive stats on the acc sampling schedules we have so far
summary(allStudies_noDups$acc_sampl_freq_per_axis)
summary(allStudies_noDups$acc_burst_duration_s)
summary(allStudies_noDups$n_samples_per_axis)

# Check the distribution of meanVedba per device type
table(allStudies_noDups$deviceType)
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="eobs"], breaks="FD", main="eobs")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="Ornitela"], breaks="FD", main="Ornitela")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="Technosmart"], breaks="FD", main="Technosmart")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="DailyDiary"], breaks="FD", main="Daily Diary")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="milsar"], breaks="FD", main="Milsar")
hist(allStudies_noDups$meanVedba_Gs[allStudies_noDups$deviceType=="madebytheo"], breaks="FD", main="madebytheo")
hist(allStudies_noDups$meanVedba_Gs, breaks="FD", main="All devices")

# Export a csv with the infos of the studies included so far, until the Vedba analysis
summaryStudyInfos <- allStudies_noDups[!duplicated(allStudies_noDups[,c("study.name","study.id","dataSource")]), c("study.name","study.id","dataSource")]
summaryStudyInfos_MB <- summaryStudyInfos[summaryStudyInfos$dataSource == "Movebank",]
write.csv(summaryStudyInfos, "DataFinalSummary/studiesIncludedInCOT_quickCheck_step3A.csv", row.names = F)
write.csv(summaryStudyInfos_MB, "DataFinalSummary/studiesIncludedInCOT_fromMovebank_quickCheck_step3A.csv", row.names = F)

# Check eobs tag generations:
onlyEobs <- allStudies_noDups[allStudies_noDups$deviceType=="eobs",]
table(onlyEobs$study.name[onlyEobs$tag.local.identifier<=2241])
