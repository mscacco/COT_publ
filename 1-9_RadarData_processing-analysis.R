
#_________________
# RADAR DATA ####
#_________________

# These script contains all analyses and processing done on the radar dataset
# (corresponding to the respective steps done on the GPA-ACC data in steps 1 to 9)
# The script can be run from line 195 using the file deposited on the Edmond repository (see readme for details)

library(data.table)
library(dplyr)

setwd("...")

#_________________________
## CALCULATE MBR and THEORETICAL COSTS based on range of body masses ###

## Import the species information including body mass and WFF_group
speciesInfos <- read.csv("DataAvailable/RADARdata_YuvalNir/midjuly_midseptember_Desert_Med_Species-Final.csv")
speciesInfos$matchingSpeciesName <- paste(speciesInfos$Genus, speciesInfos$species, sep=" ")
length(unique(speciesInfos$matchingSpeciesName))
unique(df$July)

# Check some variables
table(speciesInfos$Commonality)
table(speciesInfos$Majority.vote)
table(speciesInfos$Habitat) # All species are from all habitats (both Desert and Mediterranean)
# What about the Location? in the echo dataset I have Hazeva, Hula and Sdeboker. Are all species expected to be present at all locations?

## Associate body mass from Dunning 2007 (the same source as the rest of the data)
# Even though the radat dataset already contains body mass info, take the body mass information from the same source used for the gps data
# Import body mass infos
BMinfo <- read.csv("DataAvailable/BodyMassInfos_allBirdBatsSpecies_matchTaxonomy_gps+radar.csv", as.is=T)
# Bind, add body mass value, change some column names to match other datasets and save
table(speciesInfos$matchingSpeciesName %in% BMinfo$matchingSpeciesName)
speciesInfos_bm <- merge(speciesInfos, BMinfo[,c("matchingSpeciesName_gps_radar","species_english","BodyMass_value","BodyMass_source")], 
                         by.x="matchingSpeciesName", by.y="matchingSpeciesName_gps_radar", all.x=T)

## Reformat the month information which is one column per month. 
## The following dataset will instead have one entry per species-month combination
names(speciesInfos_bm) # month when the species is not present are ""
unique(speciesInfos_bm$July);unique(speciesInfos_bm$August);unique(speciesInfos_bm$September)
speciesMonth <- rbindlist(lapply(c("July","August","September"), function(m){
  df <- speciesInfos_bm[speciesInfos_bm[,m]!="", c("matchingSpeciesName",m,"Commonality","WFF.source","WFF","WFF_group","species_english","BodyMass_value","BodyMass_source")]
  names(df)[2] <- "month"
  df$month <- m
  return(df)
}))
table(speciesMonth$species,speciesMonth$month)

## Summarise body mass information per month and WFF_group

Mass_WFFgroup_month <- speciesMonth %>% 
  group_by(WFF_group, month) %>%
  summarise(n_species = length(unique(matchingSpeciesName)),
            meanMass = mean(BodyMass_value[Commonality %in% c("Common","common")]),
            medMass = median(BodyMass_value[Commonality %in% c("Common","common")]),
            sdMass = sd(BodyMass_value[Commonality %in% c("Common","common")]),
            minMass = min(BodyMass_value[Commonality %in% c("Common","common")]),
            maxMass = max(BodyMass_value[Commonality %in% c("Common","common")]))

saveRDS(Mass_WFFgroup_month, "DataFinalSummary/RADARdata_summaryBodyMassPerWFF-month.rds")


#_______________________
## CALCULATE OVERALL COT 
#_______________________

#________________________________
## Associate body mass to echoes

# Import body mass info per WFF and month
Mass_WFFgroup_month <- readRDS("DataFinalSummary/RADARdata_summaryBodyMassPerWFF-month.rds")

# Import Yuval's radar data with infos per echo
radarDf <- read.csv("DataAvailable/RADARdata_YuvalNir/Extended_Radar_data_midjuly_midseptember_no_bats.csv", na.strings = "NULL")
radarDf$month <- as.numeric(substr(radarDf$Time_UTC,4,5))
table(radarDf$Habitat) #all species in the species list are from all habitat, so no need to take this info into account
table(radarDf$Location)
table(radarDf$month)
radarDf$Time_UTC <- as.POSIXct(radarDf$Time_UTC, format="%d/%m/%Y %H:%M", tz="UTC")
# rename the months to match the summary df
radarDf <- radarDf %>%
  mutate(month = recode(month,
                        "7" = "July",
                        "8" = "August",
                        "9" = "September"))
table(radarDf$month)
table(Mass_WFFgroup_month$month)
table(radarDf$WFF_group)

# Add the body mass information to each echo by month and WFF group
radarDf <- radarDf %>%
  left_join(Mass_WFFgroup_month, by = c("WFF_group", "month"))
radarDf$meanMass_kg <- radarDf$meanMass/1000

# Add a variable identifying the "species" as WFFgroup-month combination
radarDf$WFFmonth_species <- paste0("WFF ",radarDf$WFF_group," - ",radarDf$month)

#_________________________
## Keep echoes < 20 sec

# Some checks on echo duration and its effect on speed
cor(radarDf$Echo_duration.s., radarDf$Speed.ms)
plot(radarDf$Echo_duration.s., radarDf$Speed.ms)
summary(radarDf$Echo_duration.s.)
hist(radarDf$Echo_duration.s.)
table(radarDf$Echo_duration.s.<20)
cor(radarDf$Echo_duration.s.[radarDf$Echo_duration.s.>20], radarDf$Speed.ms[radarDf$Echo_duration.s.>20])


# We keep for the model only echoes with duration < 20 sec
radarDf <- radarDf[radarDf$Echo_duration.s. < 20,]

cor(radarDf$Echo_duration.s., radarDf$Speed.ms)
plot(radarDf$Echo_duration.s., radarDf$Speed.ms)

#___________________________________
## Calculate metabolic rates per echo (based on MEAN body mass)

# 1. BMR and cost of passive

# BMR for birds (McKechnie & Wolf 2004, Physiological and Biochemical Zoology 77(3):502–521)
radarDf$BMR_W_mean <- 10^(-1.461 + 0.669 * log10(radarDf$meanMass)) #in grams

# MR passive flight (BMR * 2)
radarDf$theor_MR_pass_W_mean <- (2 * radarDf$BMR_W_mean)

# 2. MR flap and cost of flapping (need to run script 0C first)

# MR for flapping flight (predicted based on lm model on Kyle's data)
load("DataFinalSummary/flappingModel_KylesData.RData") #import lm object modBirds
log_bodyMass_Kg <- log(radarDf$meanMass/1000)
radarDf$theor_MR_flap_W_mean <- exp(predict(modBirds, newdata=data.frame(log_bodyMass_Kg)))

#__________________________________________________
## Calculate overall COT assuming constant flapping

# In this case, by assuming constant flapping obs_COT == theor_MR_flap 
# (because the proportion of flapping is equal one and there is no cost of passive flight to add)
# Simply transform the flap MR in W to obtain COT in J kg-1 m-1
radarDf$obs_tot_COT_J_KgM_constantFlap <-  radarDf$theor_MR_flap_W_mean / radarDf$meanMass_kg / radarDf$Speed.ms

#___________________________________________________________
## Calculate overall COT accounting for pause/pulse duration

# Calculate proportion of flapping (pulse vs pause) per echo
# column Pulse.pause is the ratio between the average pause length and average pulse length
# but to obtain the total time spent in pause and pulse we need to multiply the average pause length by the number of pauses
radarDf$tot_pauseDuration <- radarDf$Average_pause_length * radarDf$N.pauses
# we obtain the pulse duration by subtracting the pause duration from the echo duration
radarDf$tot_pulseDuration <- radarDf$Echo_duration.s. - radarDf$tot_pauseDuration
# and finally obtain the proportion of pulse (flap) per echo
radarDf$propPulse <- radarDf$tot_pulseDuration/radarDf$Echo_duration.s.
hist(radarDf$propPulse)

# Add the proportion of pulse and pause to the pass and flap MR calculation
obs_MR_pass_W_mean <- radarDf$theor_MR_pass_W_mean * (1 - radarDf$propPulse)
obs_MR_flap_W_mean <- radarDf$theor_MR_flap_W_mean * (radarDf$propPulse)

# Calculate overall flight cost per echo (in W) by summing the two
radarDf$tot_MR_W_mean <- obs_MR_pass_W_mean + obs_MR_flap_W_mean

# Remove missing speed needed for calculating final COT
radarDf <- radarDf[!is.na(radarDf$Speed.ms),]
# Calculate overall COT in J kg-1 m-1 (by dividing the tot MR by body mass in kg and by speed in m/s)
radarDf$obs_tot_COT_J_KgM_propPulse <- radarDf$tot_MR_W_mean / radarDf$meanMass_kg / radarDf$Speed.ms

# Compare the two overall COT (assuming constant flapping and accounting for the prop pulse)
boxplot(obs_tot_COT_J_KgM_propPulse~WFF_group, data=radarDf)
boxplot(obs_tot_COT_J_KgM_constantFlap~WFF_group, data=radarDf)


#_____________________________
## CALCULATE WIND SUPPORT ##
#_____________________________

source("Scripts/COT_publ/ERA5_functions_download_annotate_calculateUpliftProxies_2024Aug.R") # for wind calculation functions 

radarDf$ERA5_u_component_of_wind_ms <- as.numeric(radarDf$ERA5_u_component_of_wind_ms)
radarDf$ERA5_v_component_of_wind_ms <- as.numeric(radarDf$ERA5_v_component_of_wind_ms)

#Calculate derived wind variables from UV and save
radarDf$windSpeed_ms <- windSpeed(u=radarDf$ERA5_u_component_of_wind_ms, v=radarDf$ERA5_v_component_of_wind_ms)
radarDf$windDirFrom <- windDir_from_deg(u=radarDf$ERA5_u_component_of_wind_ms, v=radarDf$ERA5_v_component_of_wind_ms)
radarDf$windSupport <- windSupport(u=radarDf$ERA5_u_component_of_wind_ms, v=radarDf$ERA5_v_component_of_wind_ms, dg=radarDf$Azimuth)
radarDf$crossWind <- crossWind(u=radarDf$ERA5_u_component_of_wind_ms, v=radarDf$ERA5_v_component_of_wind_ms, dg=radarDf$Azimuth)
radarDf$airspeed <- airspeed(Vg=radarDf$Speed.ms, Ws=radarDf$windSupport, Cw=radarDf$crossWind)


saveRDS(radarDf, file="DataFinalSummary/RADARdata_finalSummaryDataset_perEcho_COTvariables_WFF-month_echoDurFilter.rds")
write.csv(radarDf, file="DataFinalSummary/RADARdata_finalSummaryDataset_perEcho_COTvariables_WFF-month_echoDurFilter.csv", row.names = F)

# This csv is deposited on the Edmond repository (see readme for details)

#___________
## PLOTTING ##
#___________

# Import the radar data
radarDf <- readRDS("DataFinalSummary/RADARdata_finalSummaryDataset_perEcho_COTvariables_WFF-month_echoDurFilter.rds")

radarDf$log_bodyMass_Kg <- log(radarDf$meanMass_kg)
radarDf$log_propPulse <- log(radarDf$propPulse)
radarDf$log_COT_flap <- log(radarDf$obs_tot_COT_J_KgM_constantFlap)
radarDf$log_COT_pause <- log(radarDf$obs_tot_COT_J_KgM_propPulse)

# Summary info
length(unique(radarDf$WFFmonth_species))
length(unique(radarDf$WFF_group))

radar_meanPerSp <- radarDf %>% group_by(WFFmonth_species) %>%
  summarise(log_COT_flap=mean(log_COT_flap),
            eCOT_flap=mean(obs_tot_COT_J_KgM_constantFlap),
            eCOT_pause=mean(obs_tot_COT_J_KgM_propPulse),
            log_COT_pause=mean(log_COT_pause),
            propFlap=mean(propPulse),
            log_probFlap=mean(log_propPulse),
            nSegments=n(),
            bodyMass_kg=unique(meanMass),
            log_bodyMass_Kg=unique(log_bodyMass_Kg))
radar_meanPerSp$WFFmonth_species_abb <- substr(gsub(" ","",radar_meanPerSp$WFFmonth_species),1,6)

range(radar_meanPerSp$bodyMass_kg)

#____________________________________________________
# Model and Plot prop. pulse/pause vs body mass for radar data
# The plot style below will match the plot style of the prop. flap. plot for the GPS data in script 8

radMod_pulse <- lm(log_probFlap~log_bodyMass_Kg, radar_meanPerSp)
summary(radMod_pulse)

pdf("Plots/finalPlots/radarData_propPulse_Feb2025.pdf", 5,5)
plot(log_probFlap~log_bodyMass_Kg, radar_meanPerSp, col="white", pch=20, 
     ylim=c(-0.7,0), xlim=c(-4.4,-3.85),
     ylab="log(avg Prop Pulse/Pause)", xlab="log(Body mass in Kg)",cex.lab=1.3, cex.axis=1.3)
# abline(radMod_pulse, col="grey70")
# abline(h=log(1), col="grey10") #assumes constant flapping
points(log_probFlap~log_bodyMass_Kg, radar_meanPerSp, pch=17, cex=1.8, col="grey70")
points(log_probFlap~log_bodyMass_Kg, radar_meanPerSp, pch=24, cex=1.8, col="grey10")
points(rep(log(1),nrow(radar_meanPerSp))~log_bodyMass_Kg, radar_meanPerSp, pch=8, lwd=1.5, cex=1.6, col="grey10")
text(log_probFlap~log_bodyMass_Kg, radar_meanPerSp, labels=radar_meanPerSp$WFFmonth_species_abb, pos=4, cex=1, offset=0.5, srt=30, font=3)
dev.off()


#_________________________________________
# Model and Plot COT ~ body mass for radar data
# The plot style below will match the plot style of the COT plot for the GPS data in script 9

library(lme4)
library(lmerTest)
library(ggplot2)
library(ggridges)
library(scales)

df <- rbind(data.frame(Y=radar_meanPerSp$log_COT_flap, X1=radar_meanPerSp$log_bodyMass_Kg, method="constantFlap", group=radar_meanPerSp$WFFmonth_species),
            data.frame(Y=radar_meanPerSp$log_COT_pause, X1=radar_meanPerSp$log_bodyMass_Kg, method="pulse-pause", group=radar_meanPerSp$WFFmonth_species))

length(unique(df$group))

# model
radMod_COT <- lmer(Y ~ X1 * method + (1|group), df)
(summ <- summary(radMod_COT))
r.squaredGLMM(radMod_COT)

# coefficients for reg lines
summ[["coefficients"]]
intFlap <- summ[["coefficients"]]["(Intercept)",1] # model intercept for flap
intPause <- intFlap + summ[["coefficients"]]["methodpulse-pause",1] # model intercept for swimmers
slopeFlap <- summ[["coefficients"]]["X1",1] # effect of body mass on runners
slopePause <- slopeFlap + summ[["coefficients"]]["X1:methodpulse-pause",1] # effect of body mass on swimmers
modelCoefficients <- cbind(intFlap,intPause,slopeFlap,slopePause)


# simple plot
plot(Y ~ X1, data=echoDf, col=c("gold2","dodgerblue3")[factor(echoDf$method, levels=c("constantFlap","pulse-pause"))])
abline(modelCoefficients[,"intFlap"], modelCoefficients[,"slopeFlap"], lty=1, lwd=6, col="gold2")
abline(modelCoefficients[,"intPause"], modelCoefficients[,"slopePause"], lty=1, lwd=6, col="dodgerblue3")


# Plot COT ~ body mass
load("DataFinalSummary/thereticalCOTs_Nielsen-Alexander.RData") #objects theoreticalCOTschmidt theoreticalCOTalex

pdf("Plots/finalPlots/soarersFlappers_COTmodel_Feb2025_theorLines_onlyRADAR1panel_regLines.pdf",8,5)
par(mfrow=c(1,1))
plot(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_flap, pch=19, 
     xlim=c(-4.4,-3.85), ylim=c(1,4), col="white",
     xlab="log(Body mass in Kg)", ylab="log(COT in J/ms)")
lines(log(COTrun_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex, 
      lty=2, lwd=2, col=alpha("gold2", 0.7))
lines(log(COTflight_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex,
      lty=2, lwd=2, col=alpha("forestgreen", 0.7))
lines(log(COTswim_J_KgM)~log(bodyMass_Kg), data=theoreticalCOTalex,
      lty=2, lwd=2, col=alpha("dodgerblue3", 0.7))
abline(modelCoefficients[,"intFlap"], modelCoefficients[,"slopeFlap"], lty=1, lwd=6, col="grey10")
abline(modelCoefficients[,"intPause"], modelCoefficients[,"slopePause"], lty=1, lwd=6, col="grey70")
#points(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_flap, pch=8, cex=1.8, lwd=1.5, col="grey10") # asterisks
points(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_flap, pch=17, cex=1.8, col="grey10") # triangle fill
points(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_flap, pch=24, cex=1.8, col="grey10") # triangle borders
points(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_pause, pch=17, cex=1.8, col="grey70")
points(radar_meanPerSp$log_bodyMass_Kg, radar_meanPerSp$log_COT_pause, pch=24, cex=1.8, col="grey10")
text(log_COT_flap~log_bodyMass_Kg, radar_meanPerSp, labels=radar_meanPerSp$WFFmonth_species_abb, pos=4, cex=.7, offset=0.5, srt=30, font=3)
text(log_COT_flap~log_bodyMass_Kg, radar_meanPerSp, labels=radar_meanPerSp$WFFmonth_species_abb, pos=4, cex=.7, offset=0.5, srt=30, font=3)
dev.off()

#_________________________________________
# Model residuals ~ energy landscape for radar data

# Calculate fitted (expected) COT values per group
df$Yfitted <- predict(radMod_COT, newdata = df[,c("X1","method","group")], type="response")
table(df$group)
table(df$method)

# Merge fitted value to each echo per group and method
sub1 <- radarDf[,c("log_COT_flap","log_bodyMass_Kg","WFFmonth_species","windSupport","echoID.x")]
names(sub1) <- c("Y","X1","group","windSupport","echoID")
sub1$method <- "constantFlap"
sub1 <- merge(sub1, df[df$method=="constantFlap",c("group","Yfitted")], by="group", all.x=T)
sub2 <- radarDf[,c("log_COT_pause","log_bodyMass_Kg","WFFmonth_species","windSupport","echoID.x")]
names(sub2) <- c("Y","X1","group","windSupport","echoID")
sub2$method <- "pulse-pause"
sub2 <- merge(sub2, df[df$method=="pulse-pause",c("group","Yfitted")], by="group", all.x=T)

# Calculate residuals per group from expected value
echoDf <- rbind(sub1, sub2)
echoDf$residuals <- echoDf$Y - echoDf$Yfitted

table(table(echoDf$echoID))
boxplot(residuals~method, data=echoDf)

# Model residuals ~ wind support
resMod <- lmer(residuals ~ windSupport * method + (1|echoID), data=echoDf)
(summ <- summary(resMod))
r.squaredGLMM(resMod)

# partial effect plot for residuals ~ wind support
windSupp_range <- seq(min(echoDf$windSupport, na.rm = TRUE),
                      max(echoDf$windSupport, na.rm = TRUE),
                      length.out = 20)
grid <- expand.grid(windSupport = windSupp_range,
                    method = unique(echoDf$method))
grid$COT_residuals <- predict(resMod, newdata = grid, re.form = NA)

# coeffs
intFlap <- summ[["coefficients"]]["(Intercept)",1] # model intercept for flap
intPause <- intFlap + summ[["coefficients"]]["methodpulse-pause",1] # model intercept for swimmers
slopeFlap <- summ[["coefficients"]]["windSupport",1] # effect of body mass on runners
slopePause <- slopeFlap + summ[["coefficients"]]["windSupport:methodpulse-pause",1] # effect of body mass on swimmers
modelCoefficients <- cbind(intFlap,intPause,slopeFlap,slopePause)

# simple plot of residual ~ wind support
ggplot(grid, aes(x = windSupport, y = COT_residuals, group=method, color = method)) +
  geom_line(size = 1.2) +
  labs(
    x = "Wind support",
    y = "Predicted residuals",
    color = "Method"
  ) +
  scale_color_manual(values=c("grey10","grey70")) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.15), panel.grid.major.y = element_line(linewidth = 0.15)) +
  theme(text = element_text(size = 14))
ggsave("Plots/finalPlots/radarData_residualsWindSupport.pdf", width=6,height=6)

# Plot distribution of residuals
species_order <- echoDf %>%
  group_by(group) %>%
  summarize(median_residual = median(residuals)) %>%
  arrange(desc(median_residual)) %>%
  pull(group)
echoDf$WFFmonth_species <- factor(echoDf$group, levels = species_order) 

my_gradient <- c("#0D0887","#5FC3A6","white","#B63779","#51127C")

ggplot(echoDf, aes(x = residuals, y = WFFmonth_species, fill = ..x.., group = WFFmonth_species)) +
  xlim(-3,3) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(colours = my_gradient,
                       values = scales::rescale(c(-3, 0, 3)),limits = c(-3, 3),
                       name = "COT residuals") +
  labs(x = "COT Residuals", y = "Species", fill = "COT Residuals") +
  facet_wrap(~method) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70")
ggsave("Plots/finalPlots/radar_residualsDensities_Feb2025_greenMag.pdf", width=7, height=6)

#ggplot(echoDf, aes(x = residuals, fill = method, color = method, group=method)) +
  #geom_density(alpha = 0.5, adjust = 1.25) + # For overlapping densities
ggplot(echoDf, aes(x = residuals, y = method, color = method, fill=method)) +
  geom_density_ridges(alpha = 0.5, scale = 1.1, rel_min_height = 0.01) + # For stacked densities
  xlim(-3, 3) +
  scale_color_manual(values = c("grey10", "grey70")) +
  scale_fill_manual(values = alpha(c("grey10", "grey70"), 0.7)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  labs(y = "Density", x = "Residuals")
ggsave("Plots/finalPlots/radar_residualsDensities.pdf", width=6, height=6)


