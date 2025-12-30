
setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/COT")
#setwd("/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/COT")


#____________________________
# CALCULATE COT FROM OUR DATA
#____________________________

# Import our data (one entry per segment classified as commuting)
allSegmDfs <- readRDS("DataFinalSummary/finalSummaryDataset_perSegment_fromFix_Feb2025_max10hours.rds") #dataset with one row per segment

#_______________
## Calculate BMR

# Calculate basal metabolic rate, separately for birds and bats
allSegmDfs$BMR_W <- NA
# Both of these formulas are also used by Guigueno 2019 Comparative Biochemistry and Physiology, Part A 235 193–201
# For BIRDS, phylogenetically correct (equations based on 126 avian species) (McKechnie & Wolf 2004, Physiological and Biochemical Zoology 77(3):502–521)
allSegmDfs$BMR_W[!allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")] <- 10^(-1.461 + 0.669 * log10(allSegmDfs$Body_mass_g[!allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")])) #in grams 
# For BATS from (Speakman, J.R., Thomas, D.W., Kunz, T.H., Fenton, M.B., 2003. Physiological ecology and energetics of bats. Bat ecologypp. 430–490.)
BMR_ml_O2_h <- exp(1.0895 + 0.744 * log(allSegmDfs$Body_mass_g[allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")])) #in grams (log E) 
allSegmDfs$BMR_W[allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")] <- (BMR_ml_O2_h * 20.1) / 3600 # (bmr * 20.1/3600, same as bmr * 0.005583
# For the two bats species, the BMR calculated according to this equation is slightly lower than calculated with the birds equation:
#     Species       BMR_mcKechnie   BMR_speakman
# "Pteropus lylei"    1.628944       1.203674
# "Eidolon helvum"    1.402954       1.019471

#______________________
## Calculate passive MR

# To calculate theoretical MR during passive flight we multiply the above BMR * 2
allSegmDfs$theor_MR_pass_W <- (2 * allSegmDfs$BMR_W)

#____________________________________________________
## Calculate flapping MR (need to run script 0C first)

# Predict theoretical cost of flapping based on the lm model on Kyle's data
load("DataFinalSummary/flappingModel_KylesData.RData") #modBirds, modBats, calculated in script 0C
allSegmDfs$log_bodyMass_Kg <- log(allSegmDfs$Body_mass_kg)
allSegmDfs$theor_MR_flap_W <- NA
allSegmDfs$theor_MR_flap_W[!allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")] <- exp(predict(modBirds, newdata=allSegmDfs[!allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum"),]))
allSegmDfs$theor_MR_flap_W[allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum")] <- exp(predict(modBats, newdata=allSegmDfs[allSegmDfs$species %in% c("Pteropus lylei","Eidolon helvum"),]))

#__________________________
## Calculate observed costs
## of passive and flapping, by multiplying the MR by the proportion of the two

allSegmDfs$obs_MR_pass_W <- allSegmDfs$theor_MR_pass_W *
  (allSegmDfs$avg_probPass)
allSegmDfs$obs_MR_flap_W <- allSegmDfs$theor_MR_flap_W *
  (allSegmDfs$avg_probFlap)

#_______________________
## Calculate overall COT

# Calculate overall flight cost per segment (in W) by summing the two
allSegmDfs$tot_MR_W <- allSegmDfs$obs_MR_flap_W + allSegmDfs$obs_MR_pass_W

# This should be always equal or lower (diff <= 0) than the expected metabolic rate assuming 
# constant flapping (prop. flapping == 1, meaning that observed cost == expected cost)
summary(allSegmDfs$tot_MR_W - allSegmDfs$theor_MR_flap_W)

# Calculate total cost per unit of body mass
allSegmDfs$tot_MR_W_kg <- allSegmDfs$tot_MR_W / allSegmDfs$Body_mass_kg

# Transform the observed COT from W to J kg-1 m-1
allSegmDfs$obs_tot_COT_J_KgM <- allSegmDfs$tot_MR_W / allSegmDfs$Body_mass_kg / allSegmDfs$avg_grSpeed_ms
# Transform also the theoretical flight MR to the same unit (assuming constant flapping)
allSegmDfs$theor_COT_guig_J_KgM <- allSegmDfs$theor_MR_flap_W / allSegmDfs$Body_mass_kg / allSegmDfs$avg_grSpeed_ms
# Calculate the theoretical COT of flight according to Alexander 2003
allSegmDfs$theor_COT_alex_J_KgM <- 3.6*(allSegmDfs$Body_mass_kg)^(-0.31)

# Show variation sorting species by body mass:
boxplot(log(obs_tot_COT_J_KgM)~log_bodyMass_Kg, data=allSegmDfs)

# Save the dataset with the COT info
saveRDS(allSegmDfs, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.rds")


#__________________________________________________
### Add wing information per species to the dataset

wing <- read.csv("DataFinalSummary/4supplementary_finalSpeciesList_wingMorphology_Feb2025_final.csv")

allSegmDfs_wing <- merge(allSegmDfs, wing[,c("species","wingArea_ellipse_cm2","wingLoading_kgm2")], by="species", all.x=T)

# Re-save the dataset with the wing information
saveRDS(allSegmDfs_wing, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.rds")
write.csv(allSegmDfs_wing, file="DataFinalSummary/finalSummaryDataset_perSegment_fromFix+COTvariables_Feb2025.csv", row.names = F)

#_______________________________________________
### Plot theoretical vs observed COT per species
# In this case we compare COT in Watts, between the MR predicted assuming constant flapping, and the MR calculated by summing MR pass and MR flap

library(scales)
# Aggregate infos per species
speciesCostDf <- cbind(aggregate(cbind(deviceType, theor_MR_flap_W/Body_mass_kg, Body_mass_kg)~species,
                                 allSegmDfs, FUN=unique),
                       aggregate(cbind(avg_probFlap, tot_MR_W/Body_mass_kg)~species, allSegmDfs, FUN=mean)[,2:3],
                       aggregate(Body_mass_kg~species, allSegmDfs, FUN=length)[,2],
                       aggregate(tot_MR_W/Body_mass_kg~species, allSegmDfs, FUN=IQR)[,2])
names(speciesCostDf) <- c("species","deviceType","theor_MR_flap_W_kg","Body_mass_kg","mean_flapProp","obs_COTflight_W_kg","n_obs","COT_IQR")

speciesCostDf$theor_MR_flap_W_kg <- as.numeric(speciesCostDf$theor_MR_flap_W_kg)
speciesCostDf[,c("theor_MR_flap_W_kg","mean_flapProp")] <- 
  round(speciesCostDf[,c("theor_MR_flap_W_kg","mean_flapProp")], 4)

# Order species by body mass
# speciesCostDf <- speciesCostDf[order(speciesCostDf$Body_mass_kg),]
# Order species by mean prop of flapping
speciesCostDf <- speciesCostDf[order(speciesCostDf$mean_flapProp, decreasing = T),]

# same species order in both datasets
speciesCostDf$species <- factor(speciesCostDf$species, levels = speciesCostDf$species)
allSegmDfs$species <-  factor(allSegmDfs$species, levels = speciesCostDf$species)
# device names
tagTypes <- aggregate(deviceType~species, data=allSegmDfs, FUN=function(x)paste(unique(x),collapse="|"))
tagTypes <- merge(tagTypes, data.frame(deviceType=sort(unique(tagTypes$deviceType)),
                                       deviceAbb=c("DD","E","Mt","Mt|E","M","O","T")),
                  by="deviceType", all.x=T)
tagTypes$species <- factor(tagTypes$species, levels = speciesCostDf$species)
tagTypes <- tagTypes[order(tagTypes$species),]

# Boxplot with species names
#tiff(file="Plots/finalPlots/boxplot_COTPerSpecies_ppt.tiff", 16,9, units="in", res=400)
pdf(file="Plots/finalPlots/boxplot_COTPerSpecies_ppt.pdf", 16,9)
{par(mar=c(11, 4.1, 2, 2)) #bottom, left, top, right
  bx <- boxplot(tot_MR_W_kg~species,2, data=allSegmDfs, 
                ylim=c(1,95), xaxt = "n", xlab="", outpch=20, outcol=alpha("darkgrey",0.7),
                ylab="")
  stripchart(theor_MR_flap_W_kg~species, data=speciesCostDf, col="red", 
             pch=8, vertical = TRUE, add = TRUE)
  axis(side = 1, at=1:length(bx$names), labels = FALSE)
  abline(v=1:length(bx$names), col=alpha("grey",0.4), lty=2)
  text(x=1:length(bx$names), y=speciesCostDf$theor_MR_flap_W_kg, 
       labels = speciesCostDf$n_obs, cex=0.7, pos=3, xpd = TRUE)
  text(x=1:length(bx$names), y=-5, labels = bx$names, srt = 45, adj=1, xpd = TRUE)
  text(x=1:length(bx$names), y=0, labels = tagTypes$deviceAbb, adj=0.5, xpd = TRUE, col="red", cex=0.8)
  title(xlab = substitute(paste(bold("Species (ordered by mean proportion of flapping)"))), line = 9, cex=1.2)
  title(ylab = substitute(paste(bold("Metabolic rate during flight (W Kg-1)"))), line=2.5, cex=1.2)
  #legend
  text(x=39, y=seq(67,87, length.out=6), labels=c("DD = Daily diary","E = Eobs","M = Milsar","Mt = Made by Theo","O = Ornitela","T = Technosmart-Axytreck"), adj=0, cex=0.75)
  points(x=c(39.1,39.1), y=c(90.75,94.5), col=c("grey","red"), pch=19, cex=0.8)
  text(x=39.5, y=c(90.75,94.5), labels=c("Observed","Expected"), adj=0, cex=0.75)}
dev.off()


