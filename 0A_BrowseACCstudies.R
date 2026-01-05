#_______________________________________________________________
## Browse studies available on/shared with Movebank account ####
#_______________________________________________________________

library(move)
creds <- movebankLogin()
studies <- getMovebank("study", creds)
# Sensor type IDs containing acceleration
grep("Acceleration", unique(studies$sensor_type_ids), value=T)
# Subset only those having ACC
acc_studies <- studies[grep("Acceleration", studies$sensor_type_ids),]
acc_studies$taxon_ids <- as.character(acc_studies$taxon_ids)

terrMammal <- c("Zalophus", "Vulpes", "Canis", "Wallabia", "Ziphius", "Ursus", "Lynx", "Urocyon", 
                "Homo", "Orycteropus", "Ceratotherium", "Diceros", "Capreolus", "Papio","Tremarctos","Tapirus", "Taurotragus",
                "Tapirella", "Tayassu", "Tamandua", "Syncerus", "Sus", "Suricata", "Rusa", "Spermophilus",
                "Sciurus", "Sarcophilus", "Saiga", "Saguinus", "Rhinopithecus", "Rattus", "Rangifer", "Puma","Pekania",
                "Procyon", "Lepus", "Procapra", "Prionailurus", "Priodontes", "Potos", "Pongo", "Martes", "Acinonyx",
                "Perameles", "Pecari", "Paradoxurus", "Papio", "Viverra", "Pan", "Panthera", "Chlorocebus",
                "Acinonyx", "Eira","Tragelaphus", "Loxodonta", "Equus", "Antidorcas", "Hyaena", "Lycaon", "Lycalopex",
                "Ozotoceros", "Ovis", "Orycteropus", "Odocoileus", "Ochotona", "Nycticebus", "Neovison",
                "Nasua", "Cebus", "Ateles", "Myotis", "Macrotis", "Macropus", "Macaca", "Loris", "Lophocebus",
                "Cercopithecus", "Leopardus", "Tayassu", "Lemur", "Kobus", "Giraffa", "Alcelaphus", "Connochaetes",
                "Tragelaphus", "Oryx", "Isoodon", "Ichneumia", "Hystrix", "Hydrochaeris", "Hydrochoerus", "Propithecus", 
                "Felis", "Bos", "Hapalemur", "Varecia", "Gulo", "Gorilla", "Genetta", "Oryctolagus", "Felidae", "Erinaceus",
                "Bison", "Equidae", "Elephas", "Didelphis", "Daubentonia", "Dasyprocta", "Dama", "Cuon", "Crocuta",
                "Chrysocyon", "Cerdocyon", "Chironectes", "Cervus", "Capra", "Rupicapra", "Cercocebus", "Caracal", "Canidae",
                "Budorcas", "Bubalus", "Brachylagus", "Bovidae", "Manis", "Meles", "Bassaricyon", "Antilocapra", "Ammotragus", "Alces")

waterMammal <- c("Zalophus", "Ziphius", "Tursiops", "Trichechus", "Steno", "Stenella", "Pusa", "Halichoerus",
                 "Pseudorca", "Physeter", "Phoca", "Erignathus", "Peponocephala", "Mirounga", "Grampus",
                 "Orcinus", "Arctocephalus", "Ommatophoca", "Leptonychotes", "Odobenus", "Mesoplodon",
                 "Globicephala", "Mesoplodon", "Eubalaena", "Erignathus", "Delphinus", "Cystophora", 
                 "Cetacea", "Balaenoptera", "Megaptera", "Balaena")

fish_condr <- c("Thunnus", "Sphyrna", "Spaniblennius", "Siganus", "Sciaenops", "Sphyrnidae", "Prionace", 
            "Carcharhinus", "Galeocerdo", "Isurus", "Rhincodon", "Carcharodon", "Pisodonophis", "Monodon", 
            "Mola", "Megaptera", "Manta", "Isurus", "Gadus", "Carcharias", "Anguilla", "Acipenser")

insect <- c("Exaerete frontalis", "Acherontiini")

amph <- c("Bufo", "Amphibia")

reptiles <- c("Varanus", "Tiliqua", "Testudo", "Testudinidae", "Reptilia", "Geochelone", "Chelonoidis",
              "Ophiophagus", "Mauremys", "Lepidochelys", "Kinosternon", "Indotestudo", "Gopherus", 
              "Chelonia", "Eretmochelys", "Dermochelys", "Caretta", "Bungarus", "Boiga", "Alligator", 
              "Aldabrachelys", "Agkistrodon")

bats <- c("Vespertilio", "Vespadelus", "Trachops", "Tadarida", "Sturnira", "Acerodon", 
                  "Rousettus", "Rhinolophus", "Pteropus", "Pteropodidae", "Pipistrellus",
                  "Phyllostomus", "Nyctalus", "Noctilio", "Hypsignathus", "Eidolon", "Carollia")


acc_studies_flying <- acc_studies[grep(acc_studies$taxon_ids, pattern=paste(c(terrMammal,waterMammal,fish_condr,insect,amph,reptiles),collapse="|"), invert=T),]


#____________________________
### Summary of studies #####

setwd("...")

creds <- movebankLogin() #TeamWikelski & martina.scacco

# Exploratory analysis on available data of birds:
table(acc_studies_flying$i_am_owner)
table(acc_studies_flying$i_am_collaborator)
table(acc_studies_flying$i_have_download_access)
table(acc_studies_flying$taxon_ids)

# Remove test studies
acc_studies_flying <- acc_studies_flying[!acc_studies_flying$name %in% grep("Test|test",acc_studies_flying$name, value = T),]
acc_studies_flying <- acc_studies_flying[acc_studies_flying$is_test=="false",]

# Calculate deploy duration in birds
acc_studies_flying$timestamp_first_deployed_location <- as.POSIXct(acc_studies_flying$timestamp_first_deployed_location, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
acc_studies_flying$timestamp_last_deployed_location <- as.POSIXct(acc_studies_flying$timestamp_last_deployed_location, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
acc_studies_flying$deployment_duration_days <- round(as.numeric(difftime(acc_studies_flying$timestamp_last_deployed_location, 
                                                                        acc_studies_flying$timestamp_first_deployed_location, units = "days")),1)

acc_studies_flying <- acc_studies_flying[,c("name","id",
                                          "deployment_duration_days","timestamp_first_deployed_location","timestamp_last_deployed_location","number_of_deployed_locations",
                                          "main_location_lat","main_location_long","sensor_type_ids",
                                          "i_have_download_access","i_am_owner","study_permission","principal_investigator_name",#"principal_investigator_email",
                                          "taxon_ids","license_type","license_terms")]
colnames(acc_studies_flying)[1:2] <- c("study_name","study_id")

write.csv(acc_studies_flying, file="./DataAvailable/GpsAcc_MovebankDatasets/ACC_StudiesList_BirdsBats_June2022.csv", row.names = F)

# Studies can contain more than 1 species (separated by comma)
table(acc_studies_flying$taxon_ids)

# Split by species (so that for instance a same study with 2 different species occupies 2 rows)
acc_flying_sp_allStudies <- do.call(rbind, lapply(1:nrow(acc_studies_flying), function(i){
  acc <- acc_studies_flying[i,]
  species <- strsplit(acc$taxon_ids, ",")[[1]]
  acc <- do.call(rbind, lapply(species, function(s){
    a <- acc
    if(length(strsplit(s, " ")[[1]]) == 2){
      a$genus <- strsplit(s, " ")[[1]][1]
      a$species <- s
      return(a)}
  }))
 acc <- acc[,c("study_name","study_id","genus","species",
                           "deployment_duration_days","timestamp_first_deployed_location","timestamp_last_deployed_location","number_of_deployed_locations",
                           "main_location_lat","main_location_long","sensor_type_ids",
                           "i_have_download_access","i_am_owner","study_permission","principal_investigator_name",#"principal_investigator_email",
                           "taxon_ids","license_type","license_terms")]
 return(acc)
}))

write.csv(acc_flying_sp_allStudies, file="./DataAvailable/GpsAcc_MovebankDatasets/ACC_StudiesList_BirdsBats_June2022_perSpecies.csv", row.names = F)


#___________________________
# Downloadable studies ####
# Separate those where we can actually download data (and thus the reference table)

table(acc_flying_sp_allStudies$i_am_owner)
table(acc_flying_sp_allStudies$i_have_download_access)

# change download access to TRUE for the studies accessed through my personal account
studies_martina <- c("Griffon vulture - Bardenas Reales", "Scavengers and wild ungulates", "Milvusmilvus_Milsar_SOI_final")
acc_flying_sp_allStudies$i_have_download_access[acc_flying_sp_allStudies$study_name %in% studies_martina]
acc_flying_sp_allStudies$i_have_download_access[acc_flying_sp_allStudies$study_name %in% studies_martina] <- "true"

acc_studies_download <- acc_flying_sp_allStudies[acc_flying_sp_allStudies$i_have_download_access=="true",]

# go to movebank and accept agreement for those with licence
table(acc_studies_download$license_type)
table(acc_studies_download$license_terms)
acc_studies_download$study_name[acc_studies_download$license_type=="CC_0" & acc_studies_download$i_am_owner=="false"]

length(unique(acc_studies_download$study_id))==nrow(acc_studies_download)
acc_studies_download <- acc_studies_download[!duplicated(acc_studies_download[,c("study_name","study_id","timestamp_first_deployed_location","timestamp_last_deployed_location")]),
                                             ]
write.csv(acc_studies_download, file="./DataAvailable/GpsAcc_MovebankDatasets/allACCstudies_BirdsBats_downloadPermissionT_new.csv", row.names = F)


