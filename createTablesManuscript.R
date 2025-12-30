
library(dplyr)

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
#setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")


# Summary table per study
#__________________________

df <- read.csv("DataFinalSummary/4supplementary_finalSummaryDescription_perStudy.csv")

head(df)

# Movebank studies
df1 <- df %>% 
  filter(data.source=="Movebank") %>%
  select(data.source,study.id, deviceType,species, 
              n_indiv,accSamplFreq, accBurstDurat_sec, 
              med_GPStimeLag_min) %>%
  arrange(deviceType)

# External studies
df2 <- df %>% 
  filter(data.source=="External") %>%
  arrange(deviceType) %>%
  mutate(study.id=1:10) %>%
  select(data.source,study.id, deviceType,species, 
         n_indiv,accSamplFreq, accBurstDurat_sec, 
         med_GPStimeLag_min) 

# Bind both
df_sub <- rbind(df1,df2)

# Abbreviate species column
df_sub$Species <- sapply(strsplit(df_sub$Species,"|", fixed=T), function(sps){
  paste(unlist(lapply(sps, function(sp){
  g <- paste0(substr(sp,1,1),".")
  s <- strsplit(sp, " ")[[1]][[2]]
  gs <- paste(g,s,collapse=" ")
  return(gs)
  })), collapse="|")
  })

# Round some columns
df_sub$accBurstDurat_sec <- sapply(sapply(sapply(sapply(strsplit(df_sub$accBurstDurat_sec,"|", fixed=T), as.numeric), round,1), sort), paste, collapse="|")
df_sub$med_GPStimeLag_min <- round(df_sub$med_GPStimeLag_min, 1)

# rename columns
names(df_sub) <- c("Data source","Study ID","Tag manufacturer","Species",
                   "N. individuals","ACC sampling frequency (Hz)","ACC burst duration (sec)","GPS median timelag (min)")

# Export
write.csv(df_sub, "Tables/Table_S1_studyList.csv", row.names = F)


# Summary table per species
#___________________________

df <- read.csv("DataFinalSummary/4supplementary_finalSummaryDescription_perSpecies.csv")
